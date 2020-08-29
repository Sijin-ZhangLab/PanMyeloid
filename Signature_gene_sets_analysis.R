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
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/ESCA.h5ad',"ESCA"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/HCC_10X.h5ad',"HCC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/LUNG_10X.h5ad',"LUNG"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/MEL.h5ad',"MEL"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/NPC.h5ad',"NPC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/OV-FTC.h5ad',"OV-FTC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/PAAD.h5ad',"PAAD"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/KIDNEY.h5ad',"KIDNEY"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/STAD_10X5.h5ad',"STAD"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/THCA.h5ad',"THCA"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/UCEC.h5ad',"UCEC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/MYE.h5ad',"MYE"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/LYM.h5ad',"LYM")
)


cal_DC_sig_score <- function(h5ad,Clusters,genes){
	## h5ad: the result of function parse_h5ad
	## Clusters : a vector of Clusters
	## genes : a vector of genes
	
	library(scales)
	genes_filt <- genes[genes %in% rownames(h5ad$expression)]
	expression_matrix <- h5ad$expression[genes_filt,rownames(h5ad$metadata[h5ad$metadata$MajorCluster %in% Clusters,])]
	scores <- as.data.frame(colMeans(as.matrix(expression_matrix)))
	colnames(scores) <- "value"
	scores$cluster <- h5ad$metadata[rownames(scores),]$MajorCluster
	average_score <- ddply(scores, .(cluster),function(x){mean(x$value)})
	colnames(average_score) <- c("cluster","sig_score")
	average_score$sig_score_scaled <- rescale(average_score$sig_score,to=c(0,5))
	average_score$type <- unlist(lapply(as.vector(average_score$cluster),function(x){unlist(strsplit(x,"_"))[2]}))
	type_list <- data.frame(type=c("pDC","cDC1","cDC2","cDC3"))
	average_score_merged <- merge(type_list,average_score,all.x = TRUE)
	average_score_merged$type <- factor(average_score_merged$type, levels=c("pDC","cDC1","cDC2","cDC3"))
	return(average_score_merged)
}


## score of DC signature
activeDC =c('FSCN1','BIRC3','LAMP3','CCL19','LAD1','MARCKS','TNFAIP2','CCR7','CCL22','MARCKSL1','EBI3','TNFRSF11B','NUB1','INSM1','RAB9A','LY75','SIAH2','POGLUT1','KDM2B','MGLL','TXN','MLLT6','KIF2A','GRSF1','FAM49A','PLEKHG1','SOCS2','RFTN1','AC009812.4','BMP2K','NAV1','IL7R','ID2','CCL17','PPP1R9B','NRP2','TUBB6','ARNTL2','UVRAG','TXNDC11','MREG','BTG1','NDE1','SPG11','IL32','ERICH1','TBC1D4','NFKB1','GCSAM','BZW1')

migratoryDC = c('GAL3ST','NUDT17','ITGB8','ADCY6','ENO2','IL15RA','SOCS2','IL15','STAP2','PHF24','ANKRD33B','INSM1','ANXA3','ARHGAP28','RNF115','ADORA2A','EXTL1','SPSB','SLC22A23','RABGAP1','GYG1','DAP','OGFR','GYG2','CCSER2','TMEM123','NET1','GPR52','SLCO5A1','FAH','CLU','PCGF5','SAMSN1','CDKN2B','BMP2K','ZC2HC1A','SERINC5','HIVEP1','CNR1','CNR2')

active_res <- data.frame()
mig_res <- data.frame()

for(file in h5ad_list){
	h5ad <- parse_h5ad(file[1])
	cluster <- as.vector(unique(h5ad$metadata$MajorCluster))
	cluster_used <- c(cluster[grep("pDC_LILRA4",cluster)], cluster[grep("cDC1_CLEC9A",cluster)], cluster[grep("cDC2_CD1C",cluster)], cluster[grep("cDC3_LAMP3",cluster)])
	if(length(grep("cDC3",cluster_used))){
		active_score <- cal_DC_sig_score(h5ad,cluster_used,activeDC)
		mig_score <- cal_DC_sig_score(h5ad,cluster_used,migratoryDC)
		active_score$source <- file[2]
		mig_score$source <- file[2]
		active_res <- rbind(active_res, active_score)
		mig_res <- rbind(mig_res, mig_score)
	}
	cat(file[2],"\n")
}

mig_res_formatted <- dcast(mig_res,source~type,value.var=c('sig_score_scaled'))
rownames(mig_res_formatted) <- mig_res_formatted$source
mig_res_formatted$source <- NULL
pheatmap_1 <- pheatmap::pheatmap(mig_res_formatted,show_rownames = T,show_colnames = T,cluster_cols=F,cluster_rows=F, fontsize=10,border_color ='white', na_col = "grey")


active_res_formatted <- dcast(active_res,source~type,value.var=c('sig_score_scaled'))
rownames(active_res_formatted) <- active_res_formatted$source
active_res_formatted$source <- NULL
pheatmap_2 <- pheatmap::pheatmap(active_res_formatted,show_rownames = T,show_colnames = T,cluster_cols=F,cluster_rows=F, fontsize=10,border_color ='white', na_col = "grey")

library(purrr)
library(gridExtra)
plot_list <- ls(pattern = "pheatmap_") %>% 
  map(~eval(as.name(.))[[4]])

grid.arrange(arrangeGrob(grobs = plot_list, ncol=2))


## score of macrophage signature of each dataset
angiogenesis <- c("CCND2","CCNE1","CD44","CXCR4","E2F3","EDN1","EZH2","FGF18","FGFR1","FYN","HEY1","ITGAV","JAG1","JAG2","MMP9","NOTCH1","PDGFA","PTK2","SPP1","STC1","TNFAIP6","TYMP","VAV2","VCAN","VEGFA")

phagocytosis <- c("MRC1", "CD163", "MERTK", "C1QB")

cal_Macro_sig_score <- function(h5ad,Clusters,genes){
	## h5ad: the result of function parse_h5ad
	## Clusters : a vector of Clusters
	## genes : a vector of genes
	
	library(scales)
	genes_filt <- genes[genes %in% rownames(h5ad$expression)]
	expression_matrix <- h5ad$expression[genes_filt,rownames(h5ad$metadata[h5ad$metadata$MajorCluster %in% Clusters,])]
	scores <- as.data.frame(colMeans(as.matrix(expression_matrix)))
	colnames(scores) <- "value"
	scores$cluster <- h5ad$metadata[rownames(scores),]$MajorCluster
	average_score <- ddply(scores, .(cluster),function(x){mean(x$value)})
	colnames(average_score) <- c("cluster","sig_score")
	average_score$sig_score_scaled <- rescale(average_score$sig_score,to=c(0,5))
	return(average_score)
}

angio_res <- data.frame()
phag_res <- data.frame()
meta_res <- data.frame()

for(file in h5ad_list){
	h5ad <- parse_h5ad(file[1])
	h5ad$expression <- as.matrix(h5ad$expression)
	cluster <- as.vector(unique(h5ad$metadata$MajorCluster))
	cluster_used <- c(cluster[grep("Macro",cluster)])
	#cluster_used <- cluster_used[grep('NLRP3',cluster_used,invert=TRUE)]
	cluster_used <- cluster_used[grep('LYVE1',cluster_used,invert=TRUE)]
	
	angio_score <- cal_Macro_sig_score(h5ad,cluster_used, angiogenesis)
	phag_score <- cal_Macro_sig_score(h5ad,cluster_used, phagocytosis)
	meta_score <- cal_Macro_sig_score(h5ad,cluster_used, metastasis)
	angio_score$source <- file[2]
	phag_score$source <- file[2]
	meta_score$source <- file[2]
	angio_res <- rbind(angio_res, angio_score)
	phag_res <- rbind(phag_res, phag_score)
	meta_res <- rbind(meta_res, meta_score)
	
	cat(file[2],"\n")
}

res <- data.frame(cluster=angio_res$cluster, source=angio_res$source, angio_score=angio_res$sig_score_scaled, phag_score=phag_res$sig_score_scaled,meta_score=meta_res$sig_score_scaled)
res$cluster <- paste0(unlist(lapply(as.vector(res$cluster),function(x){unlist(strsplit(x,"_"))[2]})), "_",unlist(lapply(as.vector(res$cluster),function(x){unlist(strsplit(x,"_"))[3]})))
datasource <- sort(unique(res$source))
res$meta_score <-NULL

datasource <- 
for(i in 1:length(datasource)){
	res_subset <- res[res$source == datasource[i],]
	res_subset$source <- NULL
	rownames(res_subset) <- res_subset$cluster
	res_subset$cluster <- NULL
	fig <- paste0("heatplot_",i)
	assign(fig, pheatmap(t(res_subset),show_rownames = F,show_colnames = T,cluster_cols=F,cluster_rows=F, fontsize=10,border_color ='white', na_col = "grey", main=datasource[i], legend=FALSE,cellwidth = 20, cellheight = 20,angle_col=45))
	cat( datasource[i],"\n")
}

library(purrr)
library(gridExtra)
plot_list <- ls(pattern = "heatplot_") %>% 
  map(~eval(as.name(.))[[4]])

grid.arrange(arrangeGrob(grobs = plot_list, ncol=5))



#### M1 and M2 signatures	
# signature from LZY 
M1 <- c('IL23','TNF','CXCL9','CXCL10','CXCL11','CD86','IL1A','IL1B','IL6','CCL5','IRF5','IRF1','CD40','IDO1','KYNU','CCR7')
M2 <- c('IL4R','CCL4','CCL13','CCL20','CCL17','CCL18','CCL22','CCL24','LYVE1','VEGFA','VEGFB','VEGFC','VEGFD','EGF','CTSA','CTSB','CTSC','CTSD','TGFB1','TGFB2','TGFB3','MMP14','MMP19','MMP9','CLEC7A','WNT7B','FASL','TNFSF12','TNFSF8','CD276','VTCN1','MSR1','FN1','IRF4')


cal_Macro_sig_score <- function(h5ad,Clusters,genes){
	## h5ad: the result of function parse_h5ad
	## Clusters : a vector of Clusters
	## genes : a vector of genes
	
	library(scales)
	genes_filt <- genes[genes %in% rownames(h5ad$expression)]
	expression_matrix <- h5ad$expression[genes_filt,rownames(h5ad$metadata[h5ad$metadata$MajorCluster %in% Clusters,])]
	scores <- as.data.frame(colMeans(as.matrix(expression_matrix)))
	colnames(scores) <- "value"
	scores$cluster <- h5ad$metadata[rownames(scores),]$MajorCluster
	average_score <- ddply(scores, .(cluster),function(x){mean(x$value)})
	colnames(average_score) <- c("cluster","sig_score")
	average_score$sig_score_scaled <- rescale(average_score$sig_score,to=c(0,5))
	return(average_score)
}

M1_res <- data.frame()
M2_res <- data.frame()

for(file in h5ad_list){
	h5ad <- parse_h5ad(file[1])
	cluster <- as.vector(unique(h5ad$metadata$MajorCluster))
	cluster_used <- c(cluster[c(grep("Macro",cluster), grep("Mono",cluster))])
	
	M1_score <- cal_Macro_sig_score(h5ad,cluster_used, M1)
	M2_score <- cal_Macro_sig_score(h5ad,cluster_used, M2)
	M1_score$source <- file[2]
	M2_score$source <- file[2]
	M1_res <- rbind(M1_res, M1_score)
	M2_res <- rbind(M2_res, M2_score)
	
	cat(file[2],"\n")
}

#set different color vectors for each interval
col1 = colorRampPalette(c("#e9e9e9", 'red'))(30) #set the order of greys


res <- data.frame(cluster=M1_res$cluster, source=M1_res$source, M1=M1_res$sig_score_scaled, M2=M2_res$sig_score_scaled)
res$cluster <- paste0(unlist(lapply(as.vector(res$cluster),function(x){unlist(strsplit(x,"_"))[2]})), "_",unlist(lapply(as.vector(res$cluster),function(x){unlist(strsplit(x,"_"))[3]})))
datasource <- sort(unique(res$source))
for(i in 1:length(datasource)){
	res_subset <- res[res$source == datasource[i],]
	rownames(res_subset) <- res_subset$cluster
	res_subset$source <- NULL
	res_subset$cluster <- NULL
	fig <- paste0("heatplot_",i)
	assign(fig, pheatmap(t(res_subset),show_rownames = T,show_colnames = T,cluster_cols=F,cluster_rows=F, fontsize=10,border_color ='white', na_col = "grey", main=datasource[i], color=col1,angle_col = 45,legend=FALSE,cellwidth = 20, cellheight = 20))
	cat( datasource[i],"\n")
}

library(purrr)
library(gridExtra)
plot_list <- ls(pattern = "heatplot_") %>% 
  map(~eval(as.name(.))[[4]])

grid.arrange(arrangeGrob(grobs = plot_list, ncol=5))

