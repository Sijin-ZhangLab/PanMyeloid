## /data2/csj/Pan_Myeloid/A20191105/SCENIC_analysis

## /data2/csj/tools/R-3.6.3/bin/R
.libPaths("/data2/csj/tools/Rlib3.6.3")
library(SCENIC)


### load gene expression
h5ad <- readRDS("/data2/csj/Pan_Myeloid/A20191105/NicheNet_analysis/cDC_for_NicheNet.rds")
cellInfo <- h5ad$metadata[h5ad$metadata$MajorCluster %in% c("cDC3-cDC2","cDC3-cDC1"),]
exprMat <- h5ad$expression[,rownames(cellInfo)]

setwd("/data2/csj/Pan_Myeloid/A20191105/SCENIC_analysis/LAMP3_cDC")
cellInfo <- data.frame(cellInfo)
cellTypeColumn <- "MajorCluster"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cbind(table(cellInfo$CellType))
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c(
                           "cDC3-cDC1"="darkorange", 
                           "cDC3-cDC2"="red3"
                           ))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

### Initialize settings
library(SCENIC)
org="hgnc" 
dbDir="/data2/csj/Pan_Myeloid/A20191105/SCENIC_analysis/cisTarget" # RcisTarget databases 
myDatasetTitle="SCENIC for LAMP3" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
# (Adjust minimum values according to your dataset)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
						   
exprMat_filtered <- exprMat[genesKept, ]

runCorrelation(exprMat_filtered, scenicOptions)

# export the data for GRNBoost
exportsForArboreto(exprMat_filtered, scenicOptions, dir = "int")
# TF list written as: int/1.1_inputTFs.txt
# Transposed expression matrix written as: int/1.1_exprMatrix_filtered_t.txt
### runGenie3(exprMat_filtered_log, scenicOptions) ## this step is time-consuming
### run GRNBoost in python
# source activate arboreto-env ## activate the conda environment we just created for GRNBoost (pandas version 0.23.0)
## using old version of pandas 
## pip show pandas
## pip uninstall pandas
## pip install pandas==0.23.0 --user
# python run_grnboost2.py
## pip install pandas --upgrade --users


### Build and score the GRN

GRNBoost_out <- read.table("int/net1_grn_output.tsv")
colnames(GRNBoost_out) <- c("TF", "Target", "weight")
saveRDS(GRNBoost_out,"int/1.4_GENIE3_linkList.Rds")
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat,skipHeatmap = TRUE)

# Regulators for cell type

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]


### pick significantly upregulated TF for each cluster
############## upregulated TF
TF_activity_DE <- function(regulonAUC,cellInfo,Cluster1){
	
	cell_name_1 <- rownames(cellInfo[cellInfo$CellType %in% Cluster1,])
	cell_name_2 <- rownames(cellInfo[!(cellInfo$CellType %in% Cluster1),])
	
	Expression_1 <- getAUC(regulonAUC)[,cell_name_1]
	Expression_2 <- getAUC(regulonAUC)[,cell_name_2]
	
	## log FC
	mean_c1 <- as.data.frame(rowMeans(as.matrix(Expression_1)))
	colnames(mean_c1) <- "mean_c1"
	mean_c2 <- as.data.frame(rowMeans(as.matrix(Expression_2)))
	colnames(mean_c2) <- "mean_c2"
	log2fc <- data.frame(log2fc = log2(mean_c1$mean_c1) - log2(mean_c2$mean_c2))
	rownames(log2fc) <- rownames(mean_c1)
	log2fc$gene <- rownames(log2fc)
	
	## wilcox test
	group.info <- data.frame(row.names = c(cell_name_1, cell_name_2))
	group.info[cell_name_1, "group"] <- "Group1"
	group.info[cell_name_2, "group"] <- "Group2"
	group.info[, "group"] <- factor(x = group.info[, "group"])
	data.use <- getAUC(regulonAUC)[, rownames(x = group.info), drop = FALSE]
	
	p_val <- sapply(
		X = 1:nrow(x = data.use),
		FUN = function(x) {
		return(wilcox.test(data.use[x, ] ~ group.info[, "group"])$p.value)
	})
	
	## BH correction
	adj_p_val <- p.adjust(p_val, method="BH")
	
	## DE table
	result <- data.frame(gene=log2fc$gene, log2FC=log2fc$log2fc, Pvalue=p_val, Adj_pval=adj_p_val, CellType = Cluster1)
	return(result)
}

library(AUCell)
TF_DE_res <- data.frame()
types <- unique(cellInfo$CellType)
for (i in types)
{
	Cluster1 <- i
	df<- TF_activity_DE(regulonAUC,cellInfo,Cluster1)
	TF_DE_res <- rbind(TF_DE_res,df)
}

## top 10 for each cluster
library(dplyr)
TF_DE_filt <- TF_DE_res[TF_DE_res$Adj_pval < 0.01,]
TF_DE_filt_top <- TF_DE_filt %>% group_by(CellType) %>% top_n(n = 10, wt = log2FC)
##TF_DE_filt_top <- TF_DE_filt %>% group_by(CellType) %>% top_n(n = 10, wt = -Adj_pval)

TF_list <- as.vector(TF_DE_filt_top$gene)

regulonAUC_filt <- regulonAUC[TF_list,]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC_filt)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
rownames(regulonActivity_byCellType_Scaled) <- gsub("_extended", "", rownames(regulonActivity_byCellType_Scaled))

## remove duplicated TF

p_TF <- pheatmap::pheatmap(t(regulonActivity_byCellType_Scaled[!duplicated(regulonActivity_byCellType_Scaled),]), fontsize_row=10,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-1, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color='white',cluster_cols=FALSE, cluster_rows=FALSE)
				   
p_TF
pdf(file="/data2/csj/Pan_Myeloid/A20191105/SCENIC_analysis/LAMP3_cDC/TF_regulon.pdf",width=6.73, height=3.44)	
p_TF
dev.off()
			   
pheatmap::pheatmap(regulonActivity_byCellType[!duplicated(regulonActivity_byCellType),], fontsize_row=10,
                   color=colorRampPalette(c("blue","white","red"))(100),
                   treeheight_row=10, treeheight_col=10, border_color='white',cluster_cols=FALSE, cluster_rows=FALSE)
				   

TF_DE_res$Significance <- ifelse(TF_DE_res$Adj_pval < 0.01 , "True", "False")				   
ggplot(TF_DE_res[TF_DE_res$CellType=='cDC3-cDC1',], aes(x = log2FC, y = -log10(Adj_pval))) +ylab("-Log10(adjusted Pvalue)")+ xlab("Cluster2  <-  Log2(Fold Change)  ->  Cluster1")+ 
                  geom_point(aes(color = Significance))+
                  scale_color_manual(values = c("grey","red")) + 
                  geom_text_repel(data = subset(TF_DE_res[TF_DE_res$CellType=='cDC3-cDC1',], Adj_pval < 0.01 & (log2FC>=0.5 | log2FC<=(-0.5)) ),
                aes(label = gene),size = 4, box.padding = unit(0.3, "lines"),point.padding = unit(0.2, "lines"))+
                        theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+theme(legend.position="right")+
                        ggtitle('Differentially Expressed Genes')+ theme(legend.text=element_text(size=16))+
                        theme(plot.title = element_text(hjust = 0.5))+theme(text = element_text(size=16),plot.title = element_text(size=16),
                                                                           axis.text.x = element_text(size=16),
                                                                           axis.text.y = element_text(size=16),
                                                                           axis.title.x = element_text(size=16),
                                                                           axis.title.y = element_text(size=16))