## /share/tools/R/R-3.6.0/bin/R
.libPaths("/data2/csj/tools/Rlib3.6")
library(Seurat)
library(dplyr)
library(monocle)

parse_h5ad <- function(adata){
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


raw <- parse_raw_h5ad("./cDC_monocle.h5ad")
meta <- raw$metadata
exp <- raw$expression

genes <- as.data.frame(rownames(raw$expression))
colnames(genes)<-"gene_short_name"
rownames(genes) <- genes$gene_short_name

## downsample cDC2
set.seed(5)
cDC2_cells <- rownames(meta[meta$MajorCluster == 'M04_cDC2_CD1C',])
downsample_cells <- as.vector(sample_n(as.data.frame(cDC2_cells), 700,replace = FALSE)$cDC2_cells)

other_cells <- rownames(meta[meta$MajorCluster != 'M04_cDC2_CD1C',])
exp_data_down <- exp[,c(other_cells,downsample_cells)]
meta_data_down <- meta[c(other_cells,downsample_cells),]

pd <- new("AnnotatedDataFrame", data = meta_data_down)	
fd <- new("AnnotatedDataFrame", data = genes)
HSMM <- newCellDataSet(exp_data_down,
                phenoData = pd,
                featureData = fd,
                lowerDetectionLimit = 0.5,
                expressionFamily = negbinomial.size())


HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)


#################################################picking genes using dpFeature
clustering_DEG_genes <- differentialGeneTest(HSMM, fullModelFormulaStr = '~MajorCluster', cores = 10)

HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:400]
HSMM <- setOrderingFilter(HSMM, ordering_genes = HSMM_ordering_genes)
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree')
HSMM <- orderCells(HSMM)

saveRDS(HSMM,"/data2/csj/Pan_Myeloid/A20191105/processed_data/cDC_monocle_dfFeature.rds")
saveRDS(HSMM_ordering_genes,"/data2/csj/Pan_Myeloid/A20191105/processed_data/genes_from_monocle.rds")

ClusterName_color_panel <- c(
  "M03_cDC1_CLEC9A" = "#6495ED", "M04_cDC2_CD1C" = "#339933", "M05_cDC3_LAMP3" = "#FF4500"
)

pdf(file="cDC_monocle.pdf",width=4.79, height=3.09)
plot_cell_trajectory(HSMM, color_by = "MajorCluster",cell_size=0.5) + theme(legend.position = "right") + scale_color_manual(name = "", values = ClusterName_color_panel)+guides(color = guide_legend(override.aes = list(size=5)))

dev.off()

HSMM <- readRDS("/data2/csj/Pan_Myeloid/A20191105/processed_data/cDC_monocle_dfFeature.rds")
pData(HSMM)$Pseudotime <- max(pData(HSMM)$Pseudotime) - pData(HSMM)$Pseudotime
plot_cell_trajectory(HSMM, color_by = "Pseudotime")

df <- pData(HSMM)

ggplot(df, aes(Pseudotime, colour = MajorCluster, fill=MajorCluster)) +
  geom_density(bw=0.5,size=1,alpha = 0.5)+theme_classic2()+ scale_fill_manual(name = "", values = ClusterName_color_panel)+scale_color_manual(name = "", values = ClusterName_color_panel)

diff_test_res <- differentialGeneTest(HSMM[HSMM_ordering_genes,],
              fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.001))

diff_test_res$gene <- rownames(diff_test_res)
sig_gene <- diff_test_res %>%  top_n(n = 100, wt = -qval)
sig_gene_names <- sig_gene$gene

sig_gene_names <- c("LAMP3","CCR7","FSCN1","CCL22",'IL7R','BIRC3',"FCER1A","CD1C",'FCGR2B','CLEC10A','PLAUR',"XCR1","CLEC9A",'DNASE1L3','IL1R2')
p <- plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                num_clusters = 2,
                cores = 1,
                show_rownames = T)
				
pdf("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate/Monocle_Gene.pdf",width=4.79,height=3)
 plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                num_clusters = 2,
                cores = 1,
                show_rownames = T)
dev.off()