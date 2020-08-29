
## load data
.libPaths("/data2/csj/tools/Rlib")
library(Seurat)
library(dplyr)
library(monocle)
options(stringsAsFactors=FALSE)
library(reticulate)

glm.predict <-
  function(train.data, train.group, downsample = FALSE, sample.cells = 0, genes.used = NA, test.data, test.group, alpha = 0.99, nfolds = 10) {
    ## Calculate the similarities of the train data and test data.
    ##
    ## Args:
    #' @train.data: A train data matrix with each cell in a column and each gene
    #' in a row.
    #' @train.group: A vector with the same length as the column of train.data.
    #' @downsample: Whether to sample cells in each cluster to the minimum cluster size.
    #' @sample.cells: Sample cells in each group of cells in train data, if 0 do not
    #' sample cells.
    #' @genes.used: Use a subset of genes in both the train and test data.
    #' @test.data: A test data matrix with each cell in a column and each gene
    #' in a row.
    #' @test.group: A vector with the same length as the column of train.data.
    #' @alpha: The elasticnet mixing parameter, with 0≤α≤1, passed to cv.glmnet.
    #' @nfolds: Number of folds, passed to cv.glmnet.
    ##
    ## Returns:
    ## The probability of each cell in the test.data to be predicted as each group.
    require(glmnet)
    require(ComplexHeatmap)
    glm.fits <- list()
    glm.predict <- list()
    if (length(genes.used) > 1) {
      train.data <- train.data[genes.used,]
      test.data <- test.data[genes.used,]
      if (length(genes.used) <= 50) {
        cat("There were less than 50 features used in the training data!\n")
      }
    }
    if (sample.cells == 0 & downsample) {
      sample.cells <- max(50, min(table(train.group)))
    }
    if (sample.cells > 0) {
      ngroup <- length(unique(train.group))
      if (ncol(train.data) >= sample.cells * ngroup) {
        cells_used <- c()
        for (groupi in sort(unique(train.group))) {
          if (length(which(train.group == groupi)) > sample.cells) {
            cells_used <-
              c(cells_used, sample(which(train.group == groupi), sample.cells))
          } else{
            cells_used <- c(cells_used, which(train.group == groupi))
          }
        }
        train.data <- train.data[, cells_used]
        train.group <- train.group[cells_used]
      }
    }
    for (groupi in sort(unique(train.group))) {
      fac <-  factor(train.group == groupi)
      glm.fits[[groupi]] <-
        cv.glmnet(x = t(train.data), fac, offset = getPopulationOffset(fac), 
                  family = 'binomial', intercept = FALSE, 
                  alpha = alpha, nfolds = nfolds, type.measure = 'class'
        )
      glm.predict[[groupi]] <-
        predict(
          object = glm.fits[[groupi]],
          newx = t(test.data),
          newoffset = rep(0, ncol(test.data)),
          s = 'lambda.min'
        )
    }
    glm.predict.df <- data.frame(do.call(cbind, glm.predict))
    colnames(glm.predict.df) <- sort(unique(train.group))
    glm.predict.df.prob <- (1 + exp(-glm.predict.df)) ** -1
    glm.cluster <-
      colnames(glm.predict.df.prob)[apply(glm.predict.df.prob, 1, which.max)]
    glm.predict.mean <-
      apply(glm.predict.df, 2, function(e)
        sapply(split(e, test.group), mean))
    glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
    heatmap <- Heatmap(
      t(glm.predict.mean.prob),
      name = 'Predicted\nSimilarity',
      column_title = 'test data',
      row_title = 'train data',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )
    return(
      list(
        test.group = test.group,
        logits = glm.predict.df,
        probability = glm.predict.df.prob,
        cluster = glm.cluster,
        heatmap = heatmap
      )
    )
  }

getPopulationOffset = function(y) {
  ## Calculate the offset value used in glm.predict.
  if (!is.factor(y))
    y = factor(y)
  if (length(levels(y)) != 2)
    stop("y must be a two-level factor")
  off = sum(y == levels(y)[2]) / length(y)
  off = log(off / (1 - off))
  return(rep(off, length(y)))
}

ad <- import("anndata", convert = FALSE)
### Lung Cancer
ada1 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/final_h5ad/2nd_data/LC_2nd.h5ad")
ada2 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/final_h5ad/LUNG_10X.h5ad")

ada1_df <- py_to_r(ada1$obs)$MajorCluster
ada2_df <- py_to_r(ada2$obs)$MajorCluster

gene1 <- rownames(py_to_r(ada1$var))
gene2 <- rownames(py_to_r(ada2$var))

ada1_exp <- t(py_to_r(ada1$raw$X))
rownames(ada1_exp) <- rownames(py_to_r(ada1$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

ada2_exp <- t(py_to_r(ada2$raw$X$toarray()))
rownames(ada2_exp) <- rownames(py_to_r(ada2$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

genes = c(intersect(gene1, gene2),"MARCO")

# genes = intersect(rownames(ada1_exp), rownames(ada2_exp))

ada3 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/data_for_manuscript/comparsion_between_different_datasets/Lung.h5ad")

# genes <- rownames(py_to_r(ada3$var))

res <- glm.predict(ada2_exp, ada2_df, downsample = TRUE, sample.cells = 0, genes.used = genes, ada1_exp, ada1_df, alpha = 0.99, nfolds = 10)

glm.predict.mean <-
      apply(res$logits, 2, function(e)
        sapply(split(e, res$test.group), mean))
glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1

library(circlize)
col_fun = colorRamp2(c(0,0.5, 1), c("#e9e9e9","white", "red"))

Heatmap(
      t(glm.predict.mean.prob),
	  col = col_fun,
      name = 'Predicted\nSimilarity',
      column_title = 'inDrop data',
      row_title = '10X data',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )


	
### Pancreatic adenocarcinoma
ad <- import("anndata", convert = FALSE)
ada1 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/final_h5ad/PACA.h5ad")
ada2 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/final_h5ad/2nd_data/PAAD_2nd.h5ad")

ada1_df <- py_to_r(ada1$obs)$MajorCluster
ada2_df <- py_to_r(ada2$obs)$MajorCluster

gene1 <- rownames(py_to_r(ada1$var))
gene2 <- rownames(py_to_r(ada2$var))

ada1_exp <- t(py_to_r(ada1$raw$X$toarray()))
rownames(ada1_exp) <- rownames(py_to_r(ada1$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

ada2_exp <- t(py_to_r(ada2$raw$X))
rownames(ada2_exp) <- rownames(py_to_r(ada2$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

genes = c(intersect(gene1, gene2),'FCGR3A')

res <- glm.predict(ada2_exp, ada2_df, downsample = TRUE, sample.cells = 0, genes.used = genes, ada1_exp, ada1_df, alpha = 0.99, nfolds = 10)

library(circlize)
col_fun = colorRamp2(c(0,0.5, 1), c("#e9e9e9","white", "red"))
glm.predict.mean <-
      apply(res$logits, 2, function(e)
        sapply(split(e, res$test.group), mean))
glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
PAAD <- Heatmap(
      t(glm.predict.mean.prob),
	  col = col_fun,
      name = 'Predicted\nSimilarity',
      column_title = '10X5 data',
      row_title = '10X3 data',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )

PAAD
	
## Renal cancer
ada1 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/final_h5ad/RC.h5ad")
ada2 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/final_h5ad/2nd_data/RC_2nd.h5ad")

ada1_df <- py_to_r(ada1$obs)$MajorCluster
ada2_df <- py_to_r(ada2$obs)$MajorCluster

gene1 <- rownames(py_to_r(ada1$var))
gene2 <- rownames(py_to_r(ada2$var))

ada1_exp <- t(py_to_r(ada1$raw$X$toarray()))
rownames(ada1_exp) <- rownames(py_to_r(ada1$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

ada2_exp <- t(py_to_r(ada2$raw$X$toarray()))
rownames(ada2_exp) <- rownames(py_to_r(ada2$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

genes = c(intersect(gene1, gene2))

ada3 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/data_for_manuscript/comparsion_between_different_datasets/Gastric.h5ad")
genes <- rownames(py_to_r(ada3$var))

res <- glm.predict(ada1_exp, ada1_df, downsample = TRUE, sample.cells = 0, genes.used = genes, ada2_exp, ada2_df, alpha = 0.99, nfolds = 10)

library(circlize)
col_fun = colorRamp2(c(0,0.5, 1), c("#e9e9e9","white", "red"))
glm.predict.mean <-
      apply(res$logits, 2, function(e)
        sapply(split(e, res$test.group), mean))
glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
RC <- Heatmap(
      t(glm.predict.mean.prob),
	  col = col_fun,
      name = 'Predicted\nSimilarity',
      column_title = '10X3 data',
      row_title = '10X5 data',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )

RC

## Gastric 
ada1 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/final_h5ad/2nd_data/STAD_2nd.h5ad")
ada2 <- ad$read_h5ad("/data2/csj/Pan_Myeloid/A20191105/final_h5ad/STAD_10X5.h5ad")

ada1_df <- py_to_r(ada1$obs)$MajorCluster
ada2_df <- py_to_r(ada2$obs)$MajorCluster

gene1 <- rownames(py_to_r(ada1$var))
gene2 <- rownames(py_to_r(ada2$var))

ada1_exp <- t(py_to_r(ada1$raw$X$toarray()))
rownames(ada1_exp) <- rownames(py_to_r(ada1$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

ada2_exp <- t(py_to_r(ada2$raw$X$toarray()))
rownames(ada2_exp) <- rownames(py_to_r(ada2$raw$var))
#colnames(ada1_exp) <- rownames(py_to_r(ada1$obs))

genes = c(intersect(gene1, gene2))


res <- glm.predict(ada2_exp, ada2_df, downsample = TRUE, sample.cells = 0, genes.used = genes, ada1_exp, ada1_df, alpha = 0.99, nfolds = 10)

library(circlize)
col_fun = colorRamp2(c(0,0.5, 1), c("#e9e9e9","white", "red"))
glm.predict.mean <-
      apply(res$logits, 2, function(e)
        sapply(split(e, res$test.group), mean))
glm.predict.mean.prob <- (1 + exp(-glm.predict.mean)) ** -1
STAD <- Heatmap(
      t(glm.predict.mean.prob),
	  col = col_fun,
      name = 'Predicted\nSimilarity',
      column_title = '10X3 data',
      row_title = '10X5 data',
      show_row_names = TRUE,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_title_gp = gpar(fontsize = 16),
      column_title_gp = gpar(fontsize = 16),
      row_names_gp = gpar(fontsize = 16),
      column_names_gp = gpar(fontsize = 16)
    )
STAD
