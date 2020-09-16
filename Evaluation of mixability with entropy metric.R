integrated_metadata <- read.csv("/data2/csj/Pan_Myeloid/A20191105/data_for_LZY/Scanorama_integration.csv", row.names = 1)
integrated_metadata$X <- c()
integrated_metadata$ClusterName <- all_cells_metadata[row.names(integrated_metadata),"ClusterName"]
cells_used <- sample(1:114047, 50000)
batch <- BatchEntropy(integrated_metadata[cells_used,c("UMAP1","UMAP2")], as.character(integrated_metadata[cells_used,"ClusterName"]), k.used = 80, dimension.used = "raw")
integrated_metadata_sample <- cbind(integrated_metadata[cells_used,], batch_entropy = batch$default[row.names(integrated_metadata[cells_used,])])
integrated_metadata_sample$ClusterName <- reorder(integrated_metadata_sample$ClusterName, integrated_metadata_sample$batch_entropy, mean)
save(integrated_metadata_sample, file = "./analyze/integrated_metadata_entropy.rda")
mixing_entropy <- 
  ggplot(integrated_metadata_sample, aes(x = ClusterName, y = batch_entropy)) +
  geom_boxplot(aes(fill = ClusterName), size = .3, alpha = .5, outlier.size = .5) +
  scale_fill_manual(values = ClusterName_color_panel) +
  labs(x = "", y = "Entropy") +
  theme_minimal() +
  coord_flip() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.position = "none")
ggsave(mixing_entropy, file = "./analyze/Figures/macrophage_mixing_entropy.pdf", width = 8, height = 8)


BatchEntropy <- function(input.data, group.id, k.used = 30, dimension.used = "tSNE") {
  ## Calculate the cell entropy according to the cell group.
  ##
  ## Args:
  #' @input.data: A matrix with each cell in a row. The column could be
  #' tSNE coordinate, expression matrix or even a dist format.
  #' @group.id: A vector with the same length as the row of input.data, or
  #' a list of vectors.
  #' @k.used: The k used to build a kNN-graph.
  #' @dimension.used: The method to reduce the dimension, tSNE by default,
  #' could also be PCA or raw.
  ##
  ## Returns:
  ## A vector with each cell's entropy.
  library(dbscan)
  library(entropy)
  if (dimension.used == "raw") {
    knn_graph <-
      kNN(x = input.data, k = k.used, sort = FALSE, search = "dist")
  }
  if (dimension.used == "tSNE") {
    cat(paste("Reduce the dimension using", dimension.used, "at", Sys.time(), "\n"))
    tSNE.coor <- Rtsne::Rtsne(input.data)
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <-
      kNN(x = tSNE.coor$Y, k = k.used, sort = FALSE, search = "dist")
  }
  if (dimension.used == "PCA") {
    cat(paste("Reduce the dimension using", dimension.used, "at", Sys.time(), "\n"))
    PCA.coor <- prcomp(input.data, center = FALSE)
    PCA.cumsd <- cumsum(PCA.coor$sdev) / sum(PCA.coor$sdev)
    nPCs.used <- which(PCA.cumsd > 0.9)[1]
    cat(paste("Buld a k-NN graph at", Sys.time(), "\n"))
    knn_graph <-
      kNN(x = PCA.coor$x[, 1:nPCs.used], k = k.used, sort = FALSE, search = "dist")
  }
  if (!is.list(group.id)) {
    group.id <- list(default = group.id)
  }
  cell_entropy <- list()
  for (i in names(group.id)) {
    knn_group <- matrix(group.id[[i]][knn_graph$id],
                        nrow = nrow(input.data),
                        byrow = FALSE)
    row.names(knn_group) <- row.names(input.data)
    colnames(knn_group) <- 1:k.used
    cat(paste("Calculate the cell entropy of", i, "at", Sys.time(), "\n"))
    cell_entropy[[i]] <- apply(knn_group, 1, function(x) {
      entropy(table(x))
    })
  }
  return(cell_entropy)
}