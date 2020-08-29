### TNF+ mast cells
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

h5ad_list <- list(c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerESCA_annotated.h5ad","ESCA"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerPACA_annotated.h5ad","PAAD"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerRC_annotated.h5ad","KIDNEY"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerTHCA_annotated_v3.h5ad","THCA"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerUCEC_annotated.h5ad","UCEC"),
c("/data2/csj/Pan_Myeloid/A20191105/processed_data/each_cancer/each_cancerOV_annotated.h5ad","OV-FTC"),
c("/data2/csj/Pan_Myeloid/processed_data/each_cancer_type/L_annotated_v2.h5ad","LYM"),
c("/data2/csj/Pan_Myeloid/processed_data/each_cancer_type/MM_annotated_v4.h5ad","MYE"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/Cell_BRCA_inDrop_Myeloid_re-annotated_revised.h5ad","BRCA"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/NM_Thienpont_myeloid_re-annotated_revised.h5ad","LUNG"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/Cell_Melanoma_MARS_Amit_re-annotated.h5ad","MEL"),
c("/data2/csj/Pan_Myeloid/published_data/processed_data/scHCC_Zhang_tenx_re-annotated2.h5ad","HCC"),
c("/data2/csj/Pan_Myeloid/scCRC/CLX_CellRanger/Myeloid_cells_re-annotated_P3_final_revised.h5ad","CRC"),
c("/data2/csj/Pan_Myeloid/scNPC/scNPC_re-annotated_v2.h5ad","NPC"),
c("/data2/csj/Pan_Myeloid/scGastric/5p_Myeloid_cells-annotated_revised.h5ad","STAD")
)

res <- data.frame()

for(file in h5ad_list){
	h5ad <- parse_h5ad(file[1])
	cluster <- as.vector(unique(h5ad$metadata$MajorCluster))
	cluster_used <- c(cluster[grep("Mast_KIT",cluster)])
	cell_used <- rownames(h5ad$metadata[h5ad$metadata$MajorCluster==cluster_used & h5ad$metadata$tissue == 'T',])
	TNF_expression_used <- as.data.frame(h5ad$expression["TNF",cell_used])
	colnames(TNF_expression_used) <- c("lgexp")
	TNF_frequency <- round(sum(TNF_expression_used$lgexp > 0)/length(TNF_expression_used$lgexp),3)
	VEGFA_expression_used <- as.data.frame(h5ad$expression["VEGFA",cell_used])
	colnames(VEGFA_expression_used) <- c("lgexp")
	VEGFA_frequency <- round(sum(VEGFA_expression_used$lgexp > 0)/length(VEGFA_expression_used$lgexp),3)
	Ratio <- round(TNF_frequency/VEGFA_frequency,3)
	res <- rbind(res,c(file[2],TNF_frequency,VEGFA_frequency,Ratio))
	cat(file[2],"\n")
}

colnames(res) <- c("cancer","TNF","VEGFA","Ratio")
result <- res[res$TNF != 'NaN',]
result$Ratio <- as.numeric(result$Ratio)
result$logRatio <- log2(as.numeric(result$Ratio))
result$group <- ifelse(result$logRatio < 0, "VEGFA", "TNF")

library(Polychrome)
set.seed(723451)
fifth <- createPalette(15, c("#00ffff", "#ff00ff", "#ffff00"), M=1000)

ggdotchart(result, x = "cancer", y = "logRatio",
           color = "cancer",                                # Color by groups
           palette = as.vector(fifth), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
		   add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           dot.size = 8,                                 # Large dot size
           label = "logRatio",                        # Add mpg values as dot labels
           font.label = list(color = "black", size = 9, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
           )+theme(legend.position='none')+ylab("Ratio of TNF/VEGFA+ Mast cell")+ggtitle("")+ geom_hline(yintercept = 0, linetype = 2, color = "black")

result <- result[result$cancer %in% c('BC','ESCA','Gastric','Lung','CRC','PACA','RC','UCEC','NPC'),] ## keep cancer with more than 200  	   
ggbarplot(result, x = "cancer", y = "logRatio",
          fill = "group",           # change fill color by mpg_level
          color = "white",            # Set bar border colors to white
          palette = 'jco',            # jco journal color palett. see ?ggpar
          sort.val = "des",           # Sort the value in ascending order
          sort.by.groups = FALSE,     # Don't sort inside each group
          x.text.angle = 90,          # Rotate vertically x axis texts
          ylab = FALSE,
          xlab = TRUE,
          legend.title = ""
          )+theme(legend.position='right')+coord_flip()+ylab("Ratio of TNF/VEGFA+ Mast cell")
		  
### separated by patient
res <- data.frame()
for(file in h5ad_list){
	h5ad <- parse_h5ad(file[1])
	cluster <- as.vector(unique(h5ad$metadata$MajorCluster))
	cluster_used <- c(cluster[grep("Mast_KIT",cluster)])
	patients <- unique(as.vector(h5ad$metadata$patient))
	for(PP in patients){
		cell_used <- rownames(h5ad$metadata[h5ad$metadata$MajorCluster==cluster_used & h5ad$metadata$tissue == 'T' & h5ad$metadata$patient == PP,])
		if(length(cell_used) > 0){
			TNF_expression_used <- as.data.frame(h5ad$expression["TNF",cell_used])
			colnames(TNF_expression_used) <- c("lgexp")
			TNF_frequency <- round(sum(TNF_expression_used$lgexp > 0)/length(TNF_expression_used$lgexp),3)
			VEGFA_expression_used <- as.data.frame(h5ad$expression["VEGFA",cell_used])
			colnames(VEGFA_expression_used) <- c("lgexp")
			VEGFA_frequency <- round(sum(VEGFA_expression_used$lgexp > 0)/length(VEGFA_expression_used$lgexp),3)
			Ratio <- round((TNF_frequency + 0.001)/(VEGFA_frequency+0.001),3)
			res <- rbind(res,c(file[2],PP,length(cell_used),TNF_frequency,VEGFA_frequency,Ratio))
		}
	}
	cat(file[2],"\n")
}

colnames(res) <- c("cancer","patient","cell_number","TNF","VEGFA","Ratio")
result <- res[!(res$TNF == '0' & res$VEGFA == '0'),]
result$Ratio <- as.numeric(result$Ratio)
result$logRatio <- log2(as.numeric(result$Ratio))

c54 <- c("BRCA" = 'dodgerblue2',"ESCA"='green4',"STAD"='#E31A1C',"LUNG"='#6A3D9A',"OV-FTC"='#FF7F00',
         "CRC"='#FB9A99',"PAAD"='#CAB2D6',"KIDNEY"='khaki2',"THCA"='deeppink1',"UCEC"='blue1',      
         "HCC"='steelblue4',"NPC"='green1',"MEL"='yellow4',"MYE"='yellow3',"LYM"='forestgreen')

median_table <- ddply(result,.(cancer), function(x){median(x$logRatio)})
median_table_sorted <- median_table[order(median_table$V1),]
result$cancer <- factor(result$cancer, levels=median_table_sorted$cancer)
ggboxplot(result, x = "cancer", y = "logRatio",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means()+ theme(
        legend.position = "null",
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15)
      )+xlab("")+geom_hline(yintercept=0, linetype="dashed",color = "red", size=1)+coord_flip()

res_filt <- result[result$cell_number > 10,]
median_table <- ddply(res_filt,.(cancer), function(x){median(x$logRatio)})
median_table_sorted <- median_table[order(median_table$V1),]
res_filt$cancer <- factor(res_filt$cancer, levels=median_table_sorted$cancer)
ggboxplot(res_filt, x = "cancer", y = "logRatio",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means()+ theme(
        legend.position = "null",
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),
        text = element_text(size = 10),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.text.x = element_text(
          size = 12,
          angle = 45,
          hjust = 1
        ),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15)
      )+xlab("")
	  
