## /share/tools/R/R-3.6.0/bin/R
.libPaths("/data2/csj/tools/Rlib3.6")

setwd("/data2/csj/Pan_Myeloid/A20191105/01_revised_data/RNA_bulk")
df <- read.csv("/data2/csj/Pan_Myeloid/A20191105/data_for_manuscript/umap_for_each/umap.csv")
df$ClusterName <- unlist(lapply(strsplit(as.vector(df$MajorCluster),"_"),function(x)paste0(x[2],"_",x[3])))

current.cluster.ids <- c("BRCA","CRC","ESCA","HCC","LUNG","MEL","NPC","OV-FTC","PAAD","KIDNEY","STAD","THCA","UCEC","MYE","LYM")
new.cluster.ids <- c("BRCA","CRC","ESCA","HCC","LUNG","MEL","NPC","OV","PACA","RC","STAD","THCA","UCEC","MYE","LYM")
df$cancer_new <- plyr::mapvalues(x = df$cancer, from = current.cluster.ids, to = new.cluster.ids)
df <- df[df$tissue=='T',] ## tumor tissue
df$Cancer_patient <- paste0(df$cancer_new,'.',df$patient,'.',df$tissue)

RNA_exp <- read.table("ALL.TPM.txt", header=T)
## only keep tumor tissue
RNA_exp_T <- RNA_exp[,c("symbol",colnames(RNA_exp)[grep(".T",colnames(RNA_exp))])]
## keep patients match with our pan-cancer metadata
RNA_exp_T_filt <- RNA_exp_T[,c("symbol",colnames(RNA_exp_T)[colnames(RNA_exp_T) %in% unique(df$Cancer_patient)])]

cancer_type_to_calculate <- c("ESCA","PACA","RC","THCA","UCEC")
library(plyr)
library(dplyr)
calculate_correlation <- function(cancer){
	cat(cancer,"\n")
	cancer_meta <- df[df$cancer_new==cancer & df$tissue=='T',]
	df_table <- ddply(cancer_meta,.(Cancer_patient,MajorCluster),nrow)
	colnames(df_table) <- c("patient","cluster","number")
	total <- ddply(cancer_meta,.(Cancer_patient),nrow)
	colnames(total) <- c("patient","total")

	df_total <- merge(df_table,total)
	df_total$Proportion <- round((df_total$number/df_total$total)*100, 2)
	df_total_spread <- df_total %>% select(-number,-total) %>% tidyr::spread(cluster,Proportion,fill=0) 
	rownames(df_total_spread) <- df_total_spread$patient
	df_total_spread$patient <- NULL
	#
	# extract RNA expression
	RNA_exp_T_filt_extract <- RNA_exp_T_filt[,c("symbol",colnames(RNA_exp_T_filt)[colnames(RNA_exp_T_filt) %in% rownames(df_total_spread)])]
	patient_list <- intersect(rownames(df_total_spread),colnames(RNA_exp_T_filt_extract))
	RNA_exp_T_filt_extract <- RNA_exp_T_filt_extract[,c("symbol",patient_list)]
	rownames(RNA_exp_T_filt_extract) <- RNA_exp_T_filt_extract$symbol
	RNA_exp_T_filt_extract$symbol <- NULL
	RNA_exp_T_filt_extract_t <- t(RNA_exp_T_filt_extract)
	RNA_exp_T_filt_extract_t <- RNA_exp_T_filt_extract_t[,colMeans(RNA_exp_T_filt_extract_t) >1]## we only kept genes with average TMP > 1 
	df_total_spread <- df_total_spread[patient_list,]
	
	cor_function <- function(j){
	tmp <- ldply(1:ncol(df_total_spread),function(x){cor.test(df_total_spread[,x], as.vector(RNA_exp_T_filt_extract_t[,j]), method='pearson')$p.value})
	rownames(tmp)<- colnames(df_total_spread)
	colnames(tmp) <- colnames(RNA_exp_T_filt_extract_t)[j]
	tmp <- t(tmp)
	return(tmp)
	}
	
	res <- ldply(1:ncol(RNA_exp_T_filt_extract_t),cor_function)
	rownames(res) <- colnames(RNA_exp_T_filt_extract_t)
	res <- t(res)
	res <- cbind(cancer,res)
	res <- cbind(cluster=rownames(res), res)
	return(res)
}

cor_res <- llply(cancer_type_to_calculate, calculate_correlation)
saveRDS(cor_res,"Proportion_GeneExpression_correlation.rds")

## extrac significant P values
#cor_df <- as.data.frame(cor_res[[2]]) %>% tidyr::gather(gene,pvalue,colnames(cor_res[[2]])[3]:colnames(cor_res[[2]])[ncol(cor_res[[2]])])
cor_df <- ldply(1:length(cor_res),function(x){return(as.data.frame(cor_res[[x]]) %>% tidyr::gather(gene,pvalue,colnames(cor_res[[x]])[3]:colnames(cor_res[[x]])[ncol(cor_res[[x]])]))})

result <- na.omit(cor_df)
## correct P value in each cancer

p_adjust_fun <- function(x){
	padjust <- p.adjust(x$pvalue, method ="BH")
	res <- cbind(x,padjust=padjust)
}
result_adjusted <- ddply(result,.(cancer),p_adjust_fun)

result_adjusted_filt <- result_adjusted[result_adjusted$padjust < 0.05,] ## 

result_adjusted_filt_sort <- result_adjusted_filt %>% dplyr::arrange(padjust)

## check differnt expression in Tumor an normal tissue
RNA_exp <- read.table("ALL.TPM.txt", header=T)
sample_list <- colnames(RNA_exp)
## THCA            PHF7
DE_gene <- function(cancer,gene){
	RNA_exp_cancer <- RNA_exp[,c('symbol',sample_list[grep(paste0(cancer,".P"),sample_list)])]
	RNA_exp_cancer_gene <- as.data.frame(RNA_exp_cancer[RNA_exp_cancer$symbol == gene,])
	gene_T <- as.numeric(RNA_exp_cancer_gene[1,grep(".T",colnames(RNA_exp_cancer_gene))])
	gene_N <- as.numeric(RNA_exp_cancer_gene[1,grep(".N",colnames(RNA_exp_cancer_gene))])
	pvalue <- wilcox.test(gene_T,gene_N,alternative='greater')$p.value
	log2FC <- log2(mean(gene_T)+1) - log2(mean(gene_N)+1)
	df <- data.frame(DE_log2FC=log2FC,DE_pvalue=pvalue)
	return(df)
}

add_DE_info <- function(x){
	cancer = result_adjusted_filt_sort[x,]$cancer
	gene = result_adjusted_filt_sort[x,]$gene
	DE_res <- DE_gene(cancer,gene)
	return(cbind(result_adjusted_filt_sort[x,],DE_res))
}
result_adjusted_filt_sort_addDE <- ldply(1:nrow(result_adjusted_filt_sort),function(x){add_DE_info(x)})
saveRDS(result_adjusted_filt_sort_addDE,"Proportion_GeneExpression_correlation_significant_result.rds")

## extract DE genes
final_res <- result_adjusted_filt_sort_addDE[result_adjusted_filt_sort_addDE$DE_pvalue < 0.05 & result_adjusted_filt_sort_addDE$padjust < 0.05 & result_adjusted_filt_sort_addDE$DE_log2FC > 0.5,]
final_res$name <- lapply(as.vector(final_res$cluster),function(x){paste0(strsplit(x,"_")[[1]][2:3],collapse="_")})
final_res$label <- paste0(final_res$name,"-",final_res$gene," (",final_res$cancer,")")
final_res_plot <- final_res[,c("label","padjust","DE_pvalue")]
rownames(final_res_plot) <- final_res_plot$label
final_res_plot$label <- NULL
colnames(final_res_plot) <- c("Adjusted P for correlation","P value for DE")

library(RColorBrewer)
pdf("RNA_correlation.pdf")
pheatmap::pheatmap(final_res_plot,cluster_cols=F,cluster_rows=F, fontsize=8,border_color ='white', na_col = "grey", color = colorRampPalette(rev(brewer.pal(7, "YlOrRd")))(100))
dev.off()

### visualization
library(cowplot)
library(ggplot2)
plot_scatter_plot <- function(cancer, cluster, gene){
	cancer_meta <- df[df$cancer_new==cancer & df$tissue=='T',]
	df_table <- ddply(cancer_meta,.(Cancer_patient,MajorCluster),nrow)
	colnames(df_table) <- c("patient","cluster","number")
	total <- ddply(cancer_meta,.(Cancer_patient),nrow)
	colnames(total) <- c("patient","total")

	df_total <- merge(df_table,total)
	df_total$Proportion <- round((df_total$number/df_total$total)*100, 2)
	df_total_spread <- df_total %>% select(-number,-total) %>% tidyr::spread(cluster,Proportion,fill=0) 
	rownames(df_total_spread) <- df_total_spread$patient
	df_total_spread$patient <- NULL
	#
	# extract RNA expression
	RNA_exp_T_filt_extract <- RNA_exp_T_filt[,c("symbol",colnames(RNA_exp_T_filt)[colnames(RNA_exp_T_filt) %in% rownames(df_total_spread)])]
	patient_list <- intersect(rownames(df_total_spread),colnames(RNA_exp_T_filt_extract))
	RNA_exp_T_filt_extract <- RNA_exp_T_filt_extract[,c("symbol",patient_list)]
	rownames(RNA_exp_T_filt_extract) <- RNA_exp_T_filt_extract$symbol
	RNA_exp_T_filt_extract$symbol <- NULL
	RNA_exp_T_filt_extract_t <- t(RNA_exp_T_filt_extract)
	df_total_spread <- df_total_spread[patient_list,]
	
	data_plot <- data.frame(cluster= df_total_spread[,cluster],gene=as.vector(RNA_exp_T_filt_extract_t[,gene]))
	colnames(data_plot) <- c(cluster,gene)
	
	assign("scatter_plot", 
			ggplot(data_plot, aes(x=get(gene), y=get(cluster))) + stat_cor(size=3)+geom_smooth(method='lm',formula=y~x)+ geom_point(aes(x=get(gene), y=get(cluster)),size=2)+  theme_bw(base_size = 12)+theme_cowplot()+  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background = element_blank(),panel.background = element_blank()) +theme(legend.position="right")+ xlab(gene)+ylab(cluster)
	)
  get("scatter_plot") 
}

pdf("PYCR1_SPP1.pdf",width=4,height=4)
plot_scatter_plot("THCA","M09_Macro_SPP1","PYCR1")
dev.off()

pdf("SPOCD1_SPP1.pdf")
plot_scatter_plot("THCA","M09_Macro_SPP1","SPOCD1")
dev.off()

## plot expression in normal and tumor tissue
library(ggpubr)
plot_expression <- function(cancer,gene){
	RNA_exp_cancer <- RNA_exp[,c('symbol',sample_list[grep(paste0(cancer,".P"),sample_list)])]
	rownames(RNA_exp_cancer) <- RNA_exp_cancer$symbol
	RNA_exp_cancer$symbol <- NULL
	RNA_exp_cancer_gene <- as.data.frame(t(RNA_exp_cancer[gene,]))
	RNA_exp_cancer_gene$tissue <- unlist(lapply(rownames(RNA_exp_cancer_gene),function(x){strsplit(x,'.', fixed = TRUE)[[1]][3]}))
	colnames(RNA_exp_cancer_gene) <- c("gene","tissue")
	
	assign("boxplot_plot",ggboxplot(RNA_exp_cancer_gene, x = "tissue", y = "gene",
                color = "tissue", palette =c("#00AFBB", "#E7B800"),
                add = "jitter") + stat_compare_means(label.x = 2, method.args = list(alternative="greater")) + ggtitle(paste0("Expression of ", gene," in normal and tumor tissue"))+ theme(
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
          angle = 0,
          hjust = 1
        ),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15)
      )+xlab("")
	  )
	  get("boxplot_plot") 
}

pdf("SPOCD1.pdf")
plot_expression("THCA","SPOCD1")
dev.off()

pdf("PYCR1.pdf",width=4,height=4)
plot_expression("THCA","PYCR1")
dev.off()