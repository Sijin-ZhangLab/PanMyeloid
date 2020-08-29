library(dplyr)
options(stringsAsFactors=FALSE)
library(reticulate)


ROIE <- function(crosstab){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        div.result[i,j] <- m1[i,j] / m2[i,j]
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

### updated h5ad_list
h5ad_list <- list(c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/BRCA.h5ad',"BRCA"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/CRC.h5ad',"CRC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/ESCA.h5ad',"ESCA"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/HCC.h5ad',"HCC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/LUNG.h5ad',"LUNG"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/MEL.h5ad',"MEL"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/NPC.h5ad',"NPC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/OV-FTC.h5ad',"OV-FTC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/PAAD.h5ad',"PAAD"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/KIDNEY.h5ad',"KIDNEY"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/STAD.h5ad',"STAD"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/THCA.h5ad',"THCA"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/UCEC.h5ad',"UCEC"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/MYE.h5ad',"MYE"),
c('/data2/csj/Pan_Myeloid/A20191105/final_h5ad_confirmed_by_LZY/LYM.h5ad',"LYM")
)

### the 2nd dataset
h5ad_other_list <- list(c("/data2/csj/Pan_Myeloid/scPDAC/scPDAC-annotated.h5ad","PACA"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/Science_Renal_10x_Sam_re-annotated_revised.h5ad","RC"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/Immunity_Myeloid_cells_re-annotated_revised2.h5ad","Lung"),
c("/data2/csj/Pan_Myeloid/scGastric/3p_Myeloid_cells-annotated_revised.h5ad","Gastric"))

##### tissue distribution
## only included normal and tumor tissue
res <- data.frame()
for(i in 1:length(h5ad_list)){
	h5ad <- parse_h5ad(h5ad_list[[i]][1])
	meta <- h5ad$metadata
	meta_filt <- meta[meta$tissue %in% c("N","T"),]
	meta_filt$tissue <- factor(as.vector(meta_filt$tissue),levels=c("N","T"))
	if(length(unique(meta_filt$tissue)) ==2){
		summary <- table(meta_filt[,c('MajorCluster','tissue')])
		roe <- as.data.frame(ROIE(summary))
		roe$cancer <- h5ad_list[[i]][2]
		roe$cluster <- rownames(roe)
		rownames(roe) <-  NULL
		res <- rbind(res,roe)
	}
	cat(h5ad_list[[i]][2],"\n")
}		


### monocyte tissue proportion
df <- read.csv("/data2/csj/Pan_Myeloid/A20191105/data_for_manuscript/umap_for_each/umap.csv")		   
df$cluster <- unlist(lapply(strsplit(df$MajorCluster,"_"),function(x){x[2]}))
df <- df[df$cluster=='Mono',]

df$cluster_cancer <- paste0(unlist(lapply(strsplit(df$MajorCluster,"_"),function(x){paste0(x[2],"_",x[3])})),"_",df$cancer)
df <- df[grep("CD14CD16",df$cluster_cancer,invert=T),]
df <- df[df$tissue!='L',]

ggplot(df, aes(factor(cluster_cancer)))+ geom_bar(aes(fill = tissue), position = "fill")+ xlab("")+ylab("Proportion")+theme(legend.title=element_blank(),strip.background.x = element_blank())+ scale_fill_manual(values=c("#6495ED","#FF4500","#008000"))+theme_classic2()+theme(axis.text.x=element_text(angle=45,hjust=1))


summary <- table(df[,c('cluster_cancer','tissue')])
roe <- as.data.frame(ROIE(summary))

roe$marker <- rownames(roe)
roe$cancer <- unlist(lapply(strsplit(roe$marker,"_"),function(x){x[3]}))
ggdotchart(roe, x = "marker", y = "P",
           color = "cancer",                                # Color by groups
           #palette = as.vector(fifth), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
		   add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 6,                                 # Large dot size
		   group = "cancer",
           label = round(roe$P,2),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+ geom_hline(yintercept = 1, linetype = 2, color = "black")+theme(legend.position='none')+ylab("Ro/e")+ggtitle("")+theme(axis.text.x = element_text(angle = 45,vjust = 1))
		   
		   
		   
## LAMP3 tissue distribution
library(Polychrome)
set.seed(723451)
fifth <- createPalette(15, c("#00ffff", "#ff00ff", "#ffff00"), M=1000)

LAMP3_roe <- res[grep("LAMP3",res$cluster),]
LAMP3_roe <- LAMP3_roe[LAMP3_roe$cancer %in% c('ESCA','Lung','RC','THCA','UCEC','HCC','NPC'),] ### with more than 100
ggdotchart(LAMP3_roe, x = "cancer", y = "T",
           color = "cancer",                                # Color by groups
           palette = as.vector(fifth), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
		   add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 8,                                 # Large dot size
           label = round(LAMP3_roe$T,2),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+
  geom_hline(yintercept = 1, linetype = 2, color = "black")+theme(legend.position='none')+ylab("Ro/e")+ggtitle("Tissue distribution of LAMP3+ DC")

MAST_roe <- res[grep("Mast",res$cluster),]
MAST_roe <- MAST_roe[MAST_roe$cancer %in% c('BRCA','ESCA','STAD','LUNG','CRC','PAAD','KIDNEY','UCEC','NPC'),] ### with more than 200
ggdotchart(MAST_roe, x = "cancer", y = "T",
           color = "cancer",                                # Color by groups
           #palette = as.vector(fifth), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
		   add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 8,                                 # Large dot size
           label = round(MAST_roe$T,2),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+
  geom_hline(yintercept = 1, linetype = 2, color = "black")+theme(legend.position='none')+ylab("Ro/e")+ggtitle("Tissue distribution of Mast cell")


Macro_roe <- res[grep("Macro",res$cluster),]
Macro_roe$cluster <- unlist(lapply(strsplit(as.vector(Macro_roe$cluster),"_"),function(x)paste0(x[2],"_",x[3])))

Macro_roe$marker <- paste0(Macro_roe$cancer,"_",Macro_roe$cluster)
ClusterName_color_panel <- c(
  "Mast_KIT" = "#1688A7", "pDC_LILRA4" = "#7673AE",
  "cDC1_CLEC9A" = "#b3de69", "cDC2_CD1C" = "#D195F6", "cDC3_LAMP3" = "#7E285E",
  "Mono_CD14" = "#8197FF", "Mono_CD16" = "#0911E9", "Mono_CD14CD16" = "#1FDBFE",
  "Monolike_FCN1" = "#FF9E81",
  "Macro_PPARG" = "#EF5276", "Macro_VCAN" = "#B1E7E7",
  "Macro_CX3CR1" = "#B03C0B", "Macro_FN1" = "#F39800", "Macro_GPNMB" = "#E64B35",
  "Macro_INHBA" = "#A443B2", "Macro_IL1B" = "#FFE4B5", "Macro_NLRP3" = "#FFF56A",
  "Macro_ISG15" = "#EB2C1D", "Macro_LYVE1" = "#EF5276", 
  "Macro_C1QC" = "#FD7915", "Macro_SPP1" = "#FEC718",
  "Myeloid_MKI67" = "#E43EC1"
)

ggdotchart(Macro_roe, x = "marker", y = "T",
           color = "cluster",                                # Color by groups
          palette = ClusterName_color_panel, # Custom color palette
           sorting = "descending",                       # Sort value in descending order
		   add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 6,                                 # Large dot size
		   group = "cancer",
           label = round(Macro_roe$T,2),                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+ geom_hline(yintercept = 1, linetype = 2, color = "black")+theme(legend.position='right')+ylab("Ro/e")+ggtitle("")+theme(axis.text.x = element_text(angle = 45,vjust = 1))

