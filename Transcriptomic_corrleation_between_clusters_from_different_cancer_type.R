
##########################################################################
.libPaths("/data2/csj/tools/Rlib")
library(dplyr)
library(scales)
library(reticulate)
options(stringsAsFactors=FALSE)

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


## Major cluster
h5ad <- parse_h5ad("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate/All_Myeloid.h5ad")
h5ad$metadata$Cluster <- sapply(as.vector(h5ad$metadata$MajorCluster), function(x){unlist(strsplit(x,"_"))[2]})
current.cluster.ids <- c("Mono", "Mast", "Macro", "pDC", "cDC3", "cDC2", "cDC1","Monolike", "Myeloid")
new.cluster.ids <- c("Mo/Mq", "Mast", "Mo/Mq", "pDC", "cDC", "cDC", "cDC","Mo/Mq", "Myeloid")
h5ad$metadata$MCluster <- plyr::mapvalues(x = h5ad$metadata$Cluster, from = current.cluster.ids, to = new.cluster.ids)

h5ad$metadata$type <- paste0(h5ad$metadata$MCluster,"_",h5ad$metadata$cancer)
metadata <- h5ad$metadata
expression_data_t <- as.data.frame(t(h5ad$expression))
expression <- as.data.frame(expression_data_t)
expression$anno_merge <- h5ad$metadata[rownames(expression),]$type

other_expression <- expression[c(grep("Myeloid",expression$anno_merge, invert=TRUE)),]
avg_expression <- ddply(other_expression, .(anno_merge),.fun=function(x){x$anno_merge <- NULL; colMeans(x)})
rownames(avg_expression) <- avg_expression$anno_merge
avg_expression$anno_merge <- NULL
M <- (1- cor(t(avg_expression),method="pearson"))/2
library(ape)
par(mfcol=c(1,1))
## dd <- dist(M)
hc <- hclust(as.dist(M),method="complete")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)

########################################## only DC
h5ad <- parse_h5ad("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate/adata_DC_filt_1000.h5ad")

## filt dataset (keep one dataset for each cancer)
current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
h5ad$metadata$cancer <- plyr::mapvalues(x = h5ad$metadata$cancer, from = current.cluster.ids, to = new.cluster.ids)

current.cluster.ids <- c("Mast","cDC1","cDC2","cDC3","pDC")
new.cluster.ids <- c("Mast","cDC1","cDC2","LAMP3+ cDC","pDC")
h5ad$metadata$Cluster <- plyr::mapvalues(x = h5ad$metadata$Cluster, from = current.cluster.ids, to = new.cluster.ids)

h5ad$metadata$type <- paste0(h5ad$metadata$Cluster,"_",h5ad$metadata$cancer)
metadata <- h5ad$metadata

expression_data_t <- as.data.frame(t(h5ad$expression))
expression <- as.data.frame(expression_data_t)
expression$anno_merge <- h5ad$metadata[rownames(expression),]$type

other_expression <- expression[c(grep("cDC",expression$anno_merge)),]
avg_expression <- ddply(other_expression, .(anno_merge),.fun=function(x){x$anno_merge <- NULL; colMeans(x)})
rownames(avg_expression) <- avg_expression$anno_merge
avg_expression$anno_merge <- NULL
M <- (1- cor(t(avg_expression),method="pearson"))/2
library(ape)
par(mfcol=c(1,1))
## dd <- dist(M)
hc <- hclust(as.dist(M),method="complete")
plot.phylo(as.phylo(hc), type = "fan", cex = 0.8,no.margin = FALSE)

ClusterName_color_panel <- c(
  "cDC1" = "#26933B", "cDC2" = "#F3951B", "LAMP3+ cDC" = "#E71638"
) ## "#26933B","#F3951B","#E71638"

data_plot <- as.phylo(hc)
color_names <- unlist(lapply(strsplit(as.vector(data_plot$tip.label),"_"),function(x)paste0(x[1])))
plot.phylo(data_plot, tip.col = ClusterName_color_panel[color_names], cex = 0.8,no.margin = FALSE,edge.lty = 1, edge.width =2, label.offset =0.0015,direction='downwards',srt=180,adj=1)


###################################
## only Macro and Mono
h5ad <- parse_h5ad("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate/adata_Macro_filt_1000.h5ad")

metadata <- h5ad$metadata
metadata$MajorCluster <- as.vector(metadata$MajorCluster)

## using the final label from LZY
tmp <- read.csv("/data2/csj/Pan_Myeloid/A20191105/data_for_manuscript/umap_for_each/umap.csv")
rownames(tmp) <- tmp$index
tmp <- tmp[rownames(metadata),]
metadata[rownames(tmp),]$MajorCluster <- tmp$MajorCluster

current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
metadata$cancer <- plyr::mapvalues(x = metadata$cancer, from = current.cluster.ids, to = new.cluster.ids)

metadata$Cluster <- unlist(lapply(strsplit(as.vector(metadata$MajorCluster),"_"),function(x)x[2]))
metadata$gene <- unlist(lapply(strsplit(as.vector(metadata$MajorCluster),"_"),function(x)x[3]))
metadata$Marker <- paste0(metadata$Cluster,"_",metadata$gene,"_",metadata$cancer)
expression_data_t <- as.data.frame(t(h5ad$expression))
expression <- as.data.frame(expression_data_t)
expression$anno_merge <- metadata[rownames(expression),]$Marker

other_expression <- expression[c(grep("Mono",expression$anno_merge), grep("Macro",expression$anno_merge)),]
avg_expression <- ddply(other_expression, .(anno_merge),.fun=function(x){x$anno_merge <- NULL; colMeans(x)})
rownames(avg_expression) <- avg_expression$anno_merge
avg_expression$anno_merge <- NULL

M <- (1- cor(t(avg_expression),method="pearson"))/2
library(ape)
par(mfcol=c(1,1))
## dd <- dist(M)
hc <- hclust(as.dist(M),method="complete")
data_plot <- as.phylo(hc)
color_names <- unlist(lapply(strsplit(as.vector(data_plot$tip.label),"_"),function(x)paste0(x[1],"_",x[2])))
plot.phylo(data_plot, tip.col = ClusterName_color_panel[color_names],type = "fan", cex = 0.8,no.margin = FALSE,edge.lty = 1, edge.width =2, label.offset =0.0015)
pheatmap(avg_expression,clustering_method='complete',clustering_distance_rows='correlation',border_color='white')


