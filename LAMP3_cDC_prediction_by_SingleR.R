## /share/tools/R/R-3.6.0/bin/R
.libPaths("/data2/csj/tools/Rlib3.6")
library(Seurat)
library(dplyr)
library(monocle)
options(stringsAsFactors=FALSE)
library(reticulate)
library(SingleR)
library(scran)

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

h5ad_list <- list(c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerESCA_annotated.h5ad","ESCA"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerPACA_annotated.h5ad","PACA"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerRC_annotated.h5ad","RC"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerTHCA_annotated_v3.h5ad","THCA"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerUCEC_annotated.h5ad","UCEC"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerOV_annotated.h5ad","OV"),
c("/data2/csj/Pan_Myeloid/processed_data/each_cancer_type/L_annotated_v2.h5ad","L"),
c("/data2/csj/Pan_Myeloid/processed_data/each_cancer_type/MM_annotated_v4.h5ad","MM"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/Cell_BRCA_inDrop_Myeloid_re-annotated_revised.h5ad","BRCA"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/NM_Thienpont_myeloid_re-annotated_revised.h5ad","Lung"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/Cell_Melanoma_MARS_Amit_re-annotated.h5ad","Mela"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/scHCC_Zhang_tenx_re-annotated2.h5ad","HCC"),
c("/data2/csj/Pan_Myeloid/scCRC/CLX_CellRanger/Myeloid_cells_re-annotated_P3_final_revised.h5ad","CRC"),
c("/data2/csj/Pan_Myeloid/scNPC/scNPC_re-annotated_v2.h5ad","NPC"),
c("/data2/csj/Pan_Myeloid/scGastric/5p_Myeloid_cells-annotated_revised.h5ad","Gastric")
)

res <- data.frame()
df <- data.frame()
for(file in h5ad_list){
	h5ad <- parse_h5ad(file[1])
	cluster <- as.vector(unique(h5ad$metadata$MajorCluster))
	train_clusters = c(cluster[grep("cDC1_CLEC9A",cluster)], cluster[grep("cDC2_CD1C",cluster)]) ##c('M03_cDC1_CLEC9A','M04_cDC2_CD1C')
	test_cluster = cluster[grep("cDC3_LAMP3",cluster)]  ##c('M05_cDC3_LAMP3')
	
	if(length(test_cluster) >0){ ## skip dataset without LAMP3 DC
		trainMatrix <- h5ad$expression[,rownames(h5ad$metadata[h5ad$metadata$MajorCluster %in% train_clusters,])]
		testMatrix <- h5ad$expression[,rownames(h5ad$metadata[h5ad$metadata$MajorCluster %in% test_cluster,])]

		train_label = as.vector(h5ad$metadata[h5ad$metadata$MajorCluster %in% train_clusters,]$MajorCluster)
		out <- pairwiseTTests(trainMatrix, train_label, direction="up")
		markers <- getTopMarkers(out$statistics, out$pairs, n=100)
		trained <- trainSingleR(trainMatrix, labels=train_label, genes=markers)
		pred2b <- classifySingleR(testMatrix, trained)
		
		pred_table <- data.frame(SingleR_Label = pred2b$labels)
		rownames(pred_table) <- pred2b@rownames
		pred_table$cancer <- file[2]
		df <- rbind(df,pred_table)
		
		tmp <- as.data.frame(table(pred2b$labels))
		tmp$source <- file[2]
		res <- rbind(res,tmp)
		rm(trainMatrix)
		rm(testMatrix)
	}
	rm(h5ad)
}


res$type <- unlist(lapply(as.vector(res$Var1),function(x){unlist(strsplit(x,"_"))[2]}))
library(ggpubr)
datasource <- sort(unique(res$source))
for(i in 1:length(datasource)){
	res_subset <- res[res$source == datasource[i],]
	res_subset$proportion <- paste0(res_subset$type,"\n",round(res_subset$Freq/sum(res_subset$Freq)*100,1),"%")
	fig <- paste0("p",i)
	assign(fig,ggpie(res_subset, "Freq", label = "proportion", fill = "type", color = "white", palette = c("#00AFBB", "#E7B800"))+ggtitle(datasource[i])+ theme(legend.position="none"))
}


multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14,cols=5)


