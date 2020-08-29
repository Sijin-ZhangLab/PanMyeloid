options(stringsAsFactors = F)
df <- read.csv("/data2/csj/Pan_Myeloid/A20191105/data_for_manuscript/umap_for_each/umap.csv")
df$ClusterName <- unlist(lapply(strsplit(as.vector(df$MajorCluster),"_"),function(x)paste0(x[2],"_",x[3])))
df_filt <- df[df$tissue=='T' & df$cancer !='OV-FTC' & df$cancer !='MYE' & df$cancer !='LYM' ,] ## we only keep cancer with >=5 patients

library(plyr)
calculate_mast_and_other_correlation <- function(cancer){
	cancer_meta <- df_filt[df_filt$cancer==cancer & df_filt$tissue=='T',]
	df <- ddply(cancer_meta,.(patient,ClusterName),nrow)
	colnames(df) <- c("patient","cluster","number")
	total <- ddply(cancer_meta,.(patient),nrow)
	colnames(total) <- c("patient","total")

	df_total <- merge(df,total)
	df_total$Proportion <- round((df_total$number/df_total$total)*100, 2)
	df_total_spread <- df_total %>% dplyr::select(-number,-total) %>% tidyr::spread(cluster,Proportion,fill=0) 
	rownames(df_total_spread) <- df_total_spread$patient
	df_total_spread$patient <- NULL

	corr.test(df_total_spread, adjust='none')$r
}

cancers <- unique(df_filt$cancer)
res <- lapply(cancers,calculate_mast_and_other_correlation)
names(res) <- cancers

## plot heatmap
library(pheatmap)
for(i in 1:length(cancers)){
	
	fig <- paste0("Cor_heatplot_",i)
	assign(fig, pheatmap(res[[i]],show_rownames = T,show_colnames = T,cluster_cols=T,cluster_rows=T, fontsize=9,border_color ='white', na_col = "grey", treeheight_row =20, treeheight_col=20, main=cancers[i], color = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100)))
	cat( cancers[i],"\n")
}

library(purrr)
library(gridExtra)
plot_list <- ls(pattern = "Cor_heatplot_") %>% 
  map(~eval(as.name(.))[[4]])

setwd("/data2/csj/Pan_Myeloid/A20191105/01_revised_data/TRM_comparison/")
pdf("Myeloid_subset_correlation.pdf",width=10, height=12)
grid.arrange(arrangeGrob(grobs = plot_list, ncol=3))
dev.off()

