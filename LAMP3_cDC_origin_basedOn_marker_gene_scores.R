### use gene score to predict the origin of LAMP3 DC

library(Seurat)
library(dplyr)
library(monocle)
options(stringsAsFactors=FALSE)
library(reticulate)
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

DEGene <- function(h5ad,Cluster1,Cluster2){
	## h5ad: the result of function parse_h5ad
	## Cluster1 : a vector of Cluster
	## Cluster2 : a vector of Cluster
	
	pseudocount.use =1

	cell_name_1 <- rownames(h5ad$metadata[h5ad$metadata$MajorCluster %in% Cluster1,])
	cell_name_2 <- rownames(h5ad$metadata[h5ad$metadata$MajorCluster %in% Cluster2,])
	
	Expression_1 <- h5ad$expression[,cell_name_1]
	Expression_2 <- h5ad$expression[,cell_name_2]
	
	## log FC
	mean_c1 <- as.data.frame(rowMeans(as.matrix(Expression_1)))
	colnames(mean_c1) <- "mean_c1"
	mean_c2 <- as.data.frame(rowMeans(as.matrix(Expression_2)))
	colnames(mean_c2) <- "mean_c2"
	log2fc <- data.frame(log2fc = log2(mean_c1$mean_c1 + pseudocount.use) - log2(mean_c2$mean_c2 + pseudocount.use))
	rownames(log2fc) <- rownames(mean_c1)
	log2fc$gene <- rownames(log2fc)
	
	## wilcox test
	group.info <- data.frame(row.names = c(cell_name_1, cell_name_2))
	group.info[cell_name_1, "group"] <- "Group1"
	group.info[cell_name_2, "group"] <- "Group2"
	group.info[, "group"] <- factor(x = group.info[, "group"])
	data.use <- h5ad$expression[, rownames(x = group.info), drop = FALSE]
	
	p_val <- sapply(
		X = 1:nrow(x = data.use),
		FUN = function(x) {
		return(wilcox.test(data.use[x, ] ~ group.info[, "group"])$p.value)
	})
	
	## BH correction
	adj_p_val <- p.adjust(p_val, method="BH")
	
	## DE table
	result <- data.frame(gene=log2fc$gene, log2FC=log2fc$log2fc, Pvalue=p_val, Adj_pval=adj_p_val)
	return(result)
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

for(file in h5ad_list){
	h5ad <- parse_h5ad(file[1])
	cluster <- as.vector(unique(h5ad$metadata$MajorCluster))
	cDC1_cluster = cluster[grep("cDC1_CLEC9A",cluster)] ##c('M03_cDC1_CLEC9A')
	cDC2_cluster = cluster[grep("cDC2_CD1C",cluster)] ##c('M04_cDC2_CD1C')
	cDC3_cluster = cluster[grep("cDC3_LAMP3",cluster)]  ##c('M05_cDC3_LAMP3')
	
	h5ad$expression <- as.matrix(h5ad$expression)
	if(length(cDC3_cluster) >0){ ## skip dataset without LAMP3 DC
		cDC1_cDC2 <- DEGene(h5ad,cDC1_cluster,cDC2_cluster)
		cDC3 <- DEGene(h5ad,cDC3_cluster,c(cDC1_cluster,cDC2_cluster))
		cDC1_marker <- cDC1_cDC2[cDC1_cDC2$Adj_pval < 0.05  & cDC1_cDC2$log2FC > 0.3,]$gene
		cDC2_marker <- cDC1_cDC2[cDC1_cDC2$Adj_pval < 0.05  & cDC1_cDC2$log2FC < -0.3,]$gene
		cDC3_marker <- cDC3[cDC3$Adj_pval < 0.05  & cDC3$log2FC > 0.5,]$gene
		
		total_RNA <- as.data.frame(colSums(h5ad$expression))
		colnames(total_RNA) <- 'Total'
		cDC1_sum <- as.data.frame(colSums(h5ad$expression[cDC1_marker,]))
		colnames(cDC1_sum) <- 'DC1'
		cDC2_sum <- as.data.frame(colSums(h5ad$expression[cDC2_marker,]))
		colnames(cDC2_sum) <- 'DC2'
		cDC3_sum <- as.data.frame(colSums(h5ad$expression[cDC3_marker,]))
		colnames(cDC3_sum) <- 'DC3'
		
		df <- data.frame('Total' = total_RNA$Total, 'DC1'= cDC1_sum$DC1, 'DC2'= cDC2_sum$DC2, 'DC3'= cDC3_sum$DC3)
		rownames(df) <- rownames(total_RNA)
		df$annotation <- h5ad$metadata$MajorCluster
		## only keep cDC
		df <- df[df$annotation %in% c(cDC1_cluster,cDC2_cluster,cDC3_cluster),]
		df$DC_score <- log10(df$DC2/df$Total/length(cDC2_marker)) - log10(df$DC1/df$Total/length(cDC1_marker))
		df$DC3_score <- log10(df$DC3/df$Total)
		df$cancer <- file[2]
		#ggplot(df, aes(x = DC_score, y = DC3_score))+ylab("")+ xlab("")+ geom_point(aes(color = annotation))+ theme_classic()+ scale_color_manual(values=c("#6495ED","green","#FF4500")) + guides(color = guide_legend(override.aes = list(size=5))) + geom_vline(xintercept = 0, linetype="dashed")+ ggtitle(file[2])
		
		res <- rbind(res,df)
		
	}
	rm(h5ad)
	cat(file[2]," finished\n")
}


res$Cluster <- unlist(lapply(strsplit(as.vector(res$annotation),"_"),function(x)paste0(x[2],"_",x[3])))

ggplot(res, aes(x = DC_score, y = DC3_score))+ geom_point(aes(color = Cluster))+facet_wrap(~cancer, scales='free',nrow =2)+ylab("LAMP3+ cDC score")+ xlab("cDC2-cDC1 score")+ theme_classic() + guides(color = guide_legend(override.aes = list(size=5))) + geom_vline(xintercept = 0, linetype="dashed")+ scale_color_manual(values=c("#26933B","#F3951B","#E71638"))+theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank())


	