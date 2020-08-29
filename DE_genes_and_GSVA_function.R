## load data
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

############## different expression
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


##################################################### GSVA
runGSVA <- function(h5ad,Cluster1,Cluster2,kcdf="Gaussian",AdjPvalueCutoff=0.05){
	## h5ad: the result of function parse_h5ad
	## Cluster1 : Cluster1
	## Cluster2 : Cluster2
	## kcdf: kcdf="Gaussian" for continuous and 'Poisson for integer counts'
	
	require(GSVA)
	require(GSEABase)
	require(GSVAdata)
	require(clusterProfiler)
	data(c2BroadSets)
	library(limma)

	expression <- h5ad$expression[,rownames(h5ad$metadata[h5ad$metadata$MajorCluster %in% c(Cluster1,Cluster2),])]

	## change gene symbol to geneid
	gene_entrezid <- bitr(rownames(expression), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
	expression_filt <- expression[gene_entrezid$SYMBOL,]
	rownames(expression_filt) <- gene_entrezid$ENTREZID
	expression_filt <- as.matrix(expression_filt)
	
	res_gsva <- gsva(expression_filt, c2BroadSets, parallel.sz=10,kcdf=kcdf) 

	annotation_col = data.frame(CellType = factor(h5ad$metadata[colnames(res_gsva),]$MajorCluster))
  
	rownames(annotation_col) = colnames(res_gsva)

	## using limma to conduct DE analysis
	f <- factor(annotation_col$CellType)
	design <- model.matrix(~0+f)
	colnames(design) <- c("C1","C2")
	rownames(design) <- colnames(res_gsva)

	fit <- lmFit(res_gsva, design)
	cont.matrix=makeContrasts('C1-C2',levels = design)
	fit2=contrasts.fit(fit,cont.matrix)
	fit2=eBayes(fit2)
	gs <- topTable(fit2,adjust='BH', number=Inf, p.value=AdjPvalueCutoff)
	gs$cluster <- ifelse(gs$t > 0 , "C1", "C2")
	return(gs)
}

