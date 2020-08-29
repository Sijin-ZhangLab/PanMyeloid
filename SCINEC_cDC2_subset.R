## /data2/csj/Pan_Myeloid/A20191105/SCENIC_analysis

## /data2/csj/tools/R-3.6.3/bin/R
.libPaths("/data2/csj/tools/Rlib3.6.3")
library(SCENIC)
parse_h5ad <- function(adata){
	require(reticulate)
	ad <- import("anndata", convert = FALSE)
	ada <- ad$read_h5ad(adata)
	meta <- py_to_r(ada$obs)
	if(class(ada$raw$X)[1] == "scipy.sparse.csr.csr_matrix"){
		exp <- t(py_to_r(ada$raw$X$toarray()))
	}
	else{
		exp <- t(py_to_r(ada$raw$X))
	}
	rownames(exp) <- rownames(py_to_r(ada$raw$var))
	colnames(exp) <- rownames(meta)
	return(
		list(
		metadata = meta,
		expression = exp
		)
	)
}

parse_raw_h5ad <- function(adata){
	require(reticulate)
	ad <- import("anndata", convert = FALSE)
	ada <- ad$read_h5ad(adata)
	meta <- py_to_r(ada$obs)
	if(class(ada$X)[1] == "scipy.sparse.csr.csr_matrix"){
		exp <- t(py_to_r(ada$X$toarray()))
	}
	else{
		exp <- t(py_to_r(ada$X))
	}
	rownames(exp) <- rownames(py_to_r(ada$var))
	colnames(exp) <- rownames(meta)
	return(
		list(
		metadata = meta,
		expression = exp
		)
	)
}

### load gene expression
h5ad <- parse_raw_h5ad("/data2/csj/Pan_Myeloid/A20191105/processed_data/cDC2_V2_raw_monocle.h5ad")
h5ad_scaled <- parse_h5ad("/data2/csj/Pan_Myeloid/A20191105/processed_data/cDC2_V2.h5ad")
exprMat <- h5ad$expression
cellInfo <- h5ad_scaled$metadata

cellInfo <- data.frame(cellInfo)
cellTypeColumn <- "MajorCluster"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
cbind(table(cellInfo$CellType))
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("C01_cDC2"="forestgreen", 
                           "C02_cDC2_CD1A"="darkorange", 
                           "C03_cDC2_CXCL8"="magenta4", 
                           "C04_cDC2_FCN1"="hotpink", 
                           "C05_cDC2_ISG15"="red3", 
                           "C06_cDC2_CXCL9"="skyblue"
                           ))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

### Initialize settings
library(SCENIC)
org="hgnc" 
dbDir="/data2/csj/Pan_Myeloid/A20191105/SCENIC_analysis/cisTarget" # RcisTarget databases 
myDatasetTitle="SCENIC for cDC2" # choose a name for your analysis
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
exprMat_filtered_log <- log2(exprMat_filtered+1) 

# export the data for GRNBoost
exportsForArboreto(exprMat_filtered_log, scenicOptions, dir = "int")
# TF list written as: int/1.1_inputTFs.txt
# Transposed expression matrix written as: int/1.1_exprMatrix_filtered_t.txt
### runGenie3(exprMat_filtered_log, scenicOptions) ## this step is time-consuming
### run GRNBoost in python
# source activate arboreto-env ## activate the conda environment we just created for GRNBoost (pandas version 0.23.0)
# python run_grnboost2.py


### Build and score the GRN
exprMat_log <- log2(exprMat+1)
GRNBoost_out <- read.table("int/net1_grn_output.tsv")
colnames(GRNBoost_out) <- c("TF", "Target", "weight")
saveRDS(GRNBoost_out,"int/1.4_GENIE3_linkList.Rds")
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_log,skipHeatmap = TRUE)


# Binarize activity?
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
# saveRDS(scenicOptions, file="int/scenicOptions.Rds")
#runSCENIC_4_aucell_binarize(scenicOptions)


# Regulators for cell type

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
; regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     ; function(cells) rowMeans()(regulonAUC)[,cells]))
; regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

; pheatmap::pheatmap(regulonActivity_byCellType_Scaled, fontsize_row=3, 
                   ; color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   ; treeheight_row=10, treeheight_col=10, border_color=NA)


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
types <- levels(cellInfo$CellType)
for (i in types)
{
	Cluster1 <- i
	df<- TF_activity_DE(regulonAUC,cellInfo,Cluster1)
	TF_DE_res <- rbind(TF_DE_res,df)
}

## top 10 for each cluster
library(dplyr)
TF_DE_filt <- TF_DE_res[TF_DE_res$Adj_pval < 0.01,]
TF_DE_filt_top <- TF_DE_filt %>% group_by(CellType) %>% top_n(n = 5, wt = log2FC)

as.vector(TF_DE_filt_top$gene)
TF_list <- c("HES4_extended (12g)","ZNF467_extended (13g)","KLF3 (12g)",
"IRF4 (27g)","SOX4_extended (13g)","SPI1 (694g)","CD59 (21g)", "RFX5_extended (14g)","RUNX3 (12g)",
"NR2C1_extended (11g)","NFAT5 (15g)","PRDM1 (21g)",
"BCL6_extended (13g)","EGR2 (13g)","EGR4_extended (30g)",
"FOXO1 (10g)","SP1 (31g)","CEBPA (11g)",
"KLF13 (13g)","ZBTB7A (33g)","IRF9 (122g)",
"TBL1XR1 (18g)","JAZF1 (11g)","EGR2 (13g)",
"CEBPA (11g)","IRF7 (526g)","STAT1 (1052g)",
"STAT5A (48g)","SMARCB1_extended (13g)","CEBPA (11g)"
)
regulonAUC_filt <- regulonAUC[TF_list,]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC_filt)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
rownames(regulonActivity_byCellType_Scaled) <- gsub("_extended", "", rownames(regulonActivity_byCellType_Scaled))

## remove duplicated TF

p_TF <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled[!duplicated(regulonActivity_byCellType_Scaled),], fontsize_row=10,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color='white',cluster_cols=FALSE, cluster_rows=FALSE)
				   
p_TF
pdf(file="TF_regulon.pdf",width=5.13, height=6.76)	
p_TF
dev.off()
			   
## plot expression
TF <- rownames(regulonActivity_byCellType_Scaled[!duplicated(regulonActivity_byCellType_Scaled),])
TF <- unique(unlist(lapply(TF, function(x){unlist(strsplit(x,' '))[1]})))

plot_gene_list_violin_angle(h5ad_scaled,TF)
plot_gene_list_heatmap_independently(h5ad_scaled, TF, clusters = 'all')

## plot Langerin+ DC
regulon_list <- c('SPI1 (694g)', 'RUNX3 (12g)', 'CEBPA (11g)', 'CEBPB (590g)', 'CEBPD (343g)')
regulonAUC_filt <- regulonAUC[regulon_list,]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC_filt)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
rownames(regulonActivity_byCellType_Scaled) <- gsub("_extended", "", rownames(regulonActivity_byCellType_Scaled))

## remove duplicated TF

p_TF <- pheatmap::pheatmap(regulonActivity_byCellType_Scaled[!duplicated(regulonActivity_byCellType_Scaled),], fontsize_row=10,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color='white',cluster_cols=FALSE, cluster_rows=FALSE)
				   
p_TF
pdf(file="LC_regulon.pdf",width=5.41, height=3.59)	
p_TF
dev.off()

