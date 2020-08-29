.libPaths("/data2/csj/tools/Rlib")
library(Seurat)
library(dplyr)
library(monocle)
options(stringsAsFactors=FALSE)
library(reticulate)


parse_h5ad <- function(adata){
	require(reticulate)
	ad <- import("anndata", convert = FALSE)
	ada <- ad$read_h5ad(adata)
	meta <- py_to_r(ada$obs)
	if(class(ada$raw$X)[1] == "scipy.sparse.csr.csr_matrix" | class(ada$raw$X)[1] == "scipy.sparse.csc.csc_matrix"){
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


### updated h5ad_list
h5ad_list <- list(c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/BRCA.h5ad',"BRCA"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/CRC_10X.h5ad',"CRC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/HCC_10X.h5ad',"HCC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/LUNG_10X.h5ad',"LUNG"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/MEL.h5ad',"MEL"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/NPC.h5ad',"NPC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/STAD_10X5.h5ad',"STAD")
)

ref_h5ad <- parse_h5ad("/data2/csj/Pan_Myeloid/A20191105/processed_data/cDC2_V2_final.h5ad")
cDC2_part1 <- ref_h5ad$metadata[,c("patient","cancer","MajorCluster","tissue")]

train_set <- as.data.frame(t(ref_h5ad$expression))      #construct reference set
train_set$label <- as.vector(ref_h5ad$metadata$MajorCluster)
ref_genes <- rownames(ref_h5ad$expression)
rm(ref_h5ad)

library(scibet)
cDC2_res <- data.frame()

for(file in h5ad_list){
	test_h5ad <- parse_h5ad(file[1])
	test_meta <- test_h5ad$metadata
	test_meta <- test_meta[grep("cDC2",as.vector(test_meta$MajorCluster)),]
	test_set <- t(test_h5ad$expression[,grep("cDC2",as.vector(test_h5ad$metadata$MajorCluster))])	#construct query set
	
	genes <- intersect(ref_genes, rownames(test_h5ad$expression))
	rm(test_h5ad)
	train_set_filt <- as.data.frame(train_set[,c(genes,"label")])
	test_set <- as.data.frame(test_set[,genes])
	prd <- SciBet(train_set_filt, test_set)
	test_set$label <- prd
	test_meta$MajorCluster <- prd
	test_meta <- test_meta[,c("patient","cancer","MajorCluster","tissue")]
	cDC2_res <- rbind(cDC2_res, test_meta)
	cat(file[2],"\n")
}



### calculate proportion in tumors

df <- rbind(cDC2_res,cDC2_part1)
saveRDS(df,"/data2/csj/Pan_Myeloid/A20191105/processed_data/scibet_analysis_for_cDC2/cDC2_meta_from_scibet.rds")

#### calculate the proportion of CXCL8+ cDC2 in tumor
df_filt <- df[df$tissue=='T',]

current.cluster.ids <- c("BRCA","CRC",  "HCC","LUNG","MEL","NPC","STAD","ESCA","OV","PACA","RC","THCA","UCEC","L","MM")
new.cluster.ids <- c("BRCA","CRC",  "HCC","LUNG","MEL","NPC","STAD","ESCA","OV-FTC","PAAD","KIDNEY","THCA","UCEC","LYM","MYE")
df_filt$cancer <- plyr::mapvalues(x = df_filt$cancer, from = current.cluster.ids, to = new.cluster.ids)

## barplot
c54 <- c('#1E71A9','#EF7C21','#4AA74A','#CA2828','#8A64A8','#C19A92')
ggplot(df_filt, aes(factor(cancer)))+ geom_bar(aes(fill = MajorCluster), position = "fill")+ xlab("")+ylab("Proportion")+theme(legend.title=element_blank(),strip.background.x = element_blank())+ scale_fill_manual(values=c54)+theme_classic2() +
theme(axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ))


library(plyr)
cell_num <- ddply(df_filt, .(patient,cancer), nrow)
colnames(cell_num) <- c("patient","cancer","Total")

cell_num_cluster <- ddply(df_filt, .(patient,cancer,MajorCluster), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","cluster","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer+Total~cluster,value.var=c('number'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer","Total"))

cell_num_cluster_summary_formated_melt <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$Total >=30,]
cell_num_cluster_summary_formated_melt$Proportion <- round((cell_num_cluster_summary_formated_melt$value/cell_num_cluster_summary_formated_melt$Total)*100, 2)

# Color panel -----------------
 c54 <- c("BRCA" = 'dodgerblue2',"ESCA"='green4',"STAD"='#E31A1C',"LUNG"='#6A3D9A',"OV-FTC"='#FF7F00',
         "CRC"='#FB9A99',"PAAD"='#CAB2D6',"KIDNEY"='khaki2',"THCA"='deeppink1',"UCEC"='blue1',      
         "HCC"='steelblue4',"NPC"='green1',"MEL"='yellow4',"MYE"='yellow3',"LYM"='forestgreen')

type <- "C03_cDC2_IL1B"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

library(ggpubr)
pdf("cDC2_IL1B_pro.pdf",width=4.83, height = 3.23)
ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 70, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumors"))+ theme(
        legend.position = "null",
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15)
      )
dev.off()
	  
type <- "C04_cDC2_FCN1"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

library(ggpubr)
pdf("cDC2_FCN1_pro.pdf",width=4.83, height = 3.23)
ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 70, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumors"))+ theme(
        legend.position = "null",
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15)
      )
	  
dev.off()