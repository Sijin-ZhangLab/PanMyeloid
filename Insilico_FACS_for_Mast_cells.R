library(org.Hs.eg.db)
library(ks)
library(clusterProfiler)
options(stringsAsFactors = F)

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

file <- "/data2/csj/Pan_Myeloid/scNPC/scNPC_re-annotated_v2.h5ad"
h5ad <- parse_h5ad(file)

cell_used <- rownames(h5ad$metadata[h5ad$metadata$MajorCluster=='M01_Mast_KIT'& h5ad$metadata$tissue=='T',])
exp <- t(h5ad$expression[c('TNF','VEGFA'),cell_used])

ks_density<-kde(exp)
plot(ks_density, display="filled.contour2", cont=seq(10,200,by=1), main="NPC")

##### plot TNF and VEGFA in mast cells in all datasets

h5ad_list <- list(c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerESCA_annotated.h5ad","ESCA"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerPACA_annotated.h5ad","PAAD"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerRC_annotated.h5ad","KIDNEY"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerTHCA_annotated_v3.h5ad","THCA"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerUCEC_annotated.h5ad","UCEC"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerOV_annotated.h5ad","OV-FTC"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/Cell_BRCA_inDrop_Myeloid_re-annotated_revised.h5ad","BRCA"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/NM_Thienpont_myeloid_re-annotated_revised.h5ad","LUNG"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/scHCC_Zhang_tenx_re-annotated2.h5ad","HCC"),
c("/data2/csj/Pan_Myeloid/scCRC/CLX_CellRanger/Myeloid_cells_re-annotated_P3_final_revised.h5ad","CRC"),
c("/data2/csj/Pan_Myeloid/scNPC/scNPC_re-annotated_v2.h5ad","NPC"),
c("/data2/csj/Pan_Myeloid/scGastric/5p_Myeloid_cells-annotated_revised.h5ad","STAD")
)

par(mfrow=c(3,4))
for (file in h5ad_list){
	h5ad <- parse_h5ad(file[1])
	cell_used <- rownames(h5ad$metadata[h5ad$metadata$MajorCluster=='M01_Mast_KIT'& h5ad$metadata$tissue=='T',])
	exp <- t(h5ad$expression[c('TNF','VEGFA'),cell_used])
	ks_density<-kde(exp)
	TNF_expression_used <- as.data.frame(h5ad$expression["TNF",cell_used])
	colnames(TNF_expression_used) <- c("lgexp")
	TNF_frequency <- round(sum(TNF_expression_used$lgexp > 0)/length(TNF_expression_used$lgexp),3)
	VEGFA_expression_used <- as.data.frame(h5ad$expression["VEGFA",cell_used])
	colnames(VEGFA_expression_used) <- c("lgexp")
	VEGFA_frequency <- round(sum(VEGFA_expression_used$lgexp > 0)/length(VEGFA_expression_used$lgexp),3)
	plot(ks_density, display="filled.contour2", cont=seq(10,200,by=1),main=paste0(file[2],"\n(TNF=",TNF_frequency,"; VEGFA= ",VEGFA_frequency,")"))
}

