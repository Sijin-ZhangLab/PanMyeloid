library(plyr)
library(ggplot2)
library(scater)
rm(list=ls())

## plot cell and patient number
summary <- read.csv("/data2/csj/Pan_Myeloid/A20191105/data_summary_for_plot.csv")
p1<- ggbarplot(summary, "Cancer", "Patient",
   fill = "Cancer", color = "Cancer",
   palette = c54)+ theme(
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
      )+ xlab("")+ylab("Number of patients")
	  
p2<- ggbarplot(summary, "Cancer", "Cell",
   fill = "Cancer", color = "Cancer",
   palette = c54)+ theme(
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
      )+ xlab("")+ylab("Number of cells")
	  
multiplot(p1,p2,cols=1)
   
setwd("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate")
metadata <- read.csv("metadata.csv")

patient_summary <- metadata[,c("cancer","patient")]
patient_summary <- patient_summary[!duplicated(patient_summary),]
table(patient_summary$cancer)

metadata$sample <- paste0(metadata$patient,"-",metadata$tissue)
sample_summary <- metadata[,c("cancer","sample")]
sample_summary <- sample_summary[!duplicated(sample_summary),]
table(sample_summary$cancer)

metadata_filt <- metadata[metadata$Cluster!='Myeloid',]
metadata_filt <- metadata_filt[metadata_filt$tissue!='L',]
metadata_filt$tissue <- factor(metadata_filt$tissue, levels=c("T","N","P"))
metadata_filt$cancer <- factor(metadata_filt$cancer, levels=c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L"))

metadata_filt$Cluster <- factor(metadata_filt$Cluster, levels = c("Mast","pDC",  "cDC", "Mo/Mq"))

## cluster distribution
ggplot(metadata_filt, aes(factor(tissue)))+ geom_bar(aes(fill = Cluster), position = "fill")+ xlab("")+ylab("Proportion")+facet_wrap(~cancer,ncol=5,scales='free_x')+theme(legend.title=element_blank(),strip.background.x = element_blank())+ scale_fill_manual(values=c("#6495ED","#FFA500","#FF4500","#999999"))+theme_classic2()


## MajorCluster comparison
## only include N and T
metadata_filt_N_T <- metadata_filt[metadata_filt$tissue %in% c("N","T"),]

## Mast
metadata_filt_N_T <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c("BC","CRC","ESCA","Gastric","HCC","Lung","OV","PACA","RC","THCA","UCEC"),]
metadata_filt_N_T$tissue <- factor(metadata_filt_N_T$tissue)
cell_num <- ddply(metadata_filt_N_T, .(patient,cancer,tissue), nrow)
colnames(cell_num) <- c("patient","cancer","tissue","Total")

metadata_filt_Mast <- metadata_filt_N_T[metadata_filt_N_T$Cluster %in% c("Mast"),]
cell_num_cluster <- ddply(metadata_filt_Mast, .(patient,cancer,tissue), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","tissue","number")
cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

## paired wilcox test
library(ggpubr)
cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer~tissue,value.var=c('Proportion'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer"))

p <- ggpaired(cell_num_cluster_summary_formated_melt, x = "variable", y = "value",
          color='variable',palette = "jco", 
          line.color = "gray", line.size = 0.4,
          short.panel.labs = TRUE)+facet_wrap(~cancer, scales = "free",strip.position="top",ncol =3)+ ylab("Proportion of mast cells")+ xlab("")
p + stat_compare_means(paired = TRUE,method.args = list(alternative = "less"))
p + stat_compare_means(paired = TRUE,method.args = list(alternative = "greater"))

## cDC
metadata_filt_N_T <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c("BC","CRC","ESCA","Gastric","HCC","Lung","OV","PACA","RC","THCA","UCEC"),]
metadata_filt_N_T$tissue <- factor(metadata_filt_N_T$tissue)
cell_num <- ddply(metadata_filt_N_T, .(patient,cancer,tissue), nrow)
colnames(cell_num) <- c("patient","cancer","tissue","Total")

metadata_filt_cDC <- metadata_filt_N_T[metadata_filt_N_T$Cluster %in% c("cDC"),]
cell_num_cluster <- ddply(metadata_filt_cDC, .(patient,cancer,tissue), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","tissue","number")
cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=100,]
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

## paired wilcox test
library(ggpubr)
cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer~tissue,value.var=c('Proportion'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer"))

p <- ggpaired(cell_num_cluster_summary_formated_melt, x = "variable", y = "value",
          color='variable',palette = "jco", 
          line.color = "gray", line.size = 0.4,
          short.panel.labs = TRUE)+facet_wrap(~cancer, scales = "free",strip.position="top",ncol =3)+ ylab("Proportion of cDCs")+ xlab("")
p + stat_compare_means(paired = TRUE,method.args = list(alternative = "less"))
p + stat_compare_means(paired = TRUE,method.args = list(alternative = "greater"))



## DC subset comparison
## DC subset cluster distribution
metadata_filt_DC <- metadata_filt[metadata_filt$Cluster %in% c("cDC"),]
metadata_filt_DC$DC_type <- unlist(lapply(strsplit(metadata_filt_DC$MajorCluster,"_"),function(x){x[2]}))

ggplot(metadata_filt_DC, aes(factor(tissue)))+ geom_bar(aes(fill = DC_type), position = "fill")+ xlab("")+ylab("Proportion")+facet_wrap(~cancer,ncol=5,scales='free_x')+theme(legend.title=element_blank(),strip.background.x = element_blank())+scale_fill_brewer(palette="Set1")

metadata_filt_DC$DC_type <- factor(metadata_filt_DC$DC_type, levels=c("pDC",'cDC1','cDC2','cDC3'))
current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
metadata_filt_DC$cancer <- plyr::mapvalues(x = metadata_filt_DC$cancer, from = current.cluster.ids, to = new.cluster.ids)
ggplot(metadata_filt_DC, aes(factor(tissue)))+ geom_bar(aes(fill = DC_type), position = "fill")+ xlab("")+ylab("Proportion")+facet_wrap(~cancer,ncol=5,scales='free_x')+theme(legend.title=element_blank(),strip.background.x = element_blank())+ scale_fill_manual(values=c("#26933B","#F3951B","#E71638"))+theme_classic2()

## LAMP3 DC
metadata_filt_DC_N_T <- metadata_filt_DC[metadata_filt_DC$tissue %in% c("N","T"),]
current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
metadata_filt_DC_N_T$cancer <- plyr::mapvalues(x = metadata_filt_DC_N_T$cancer, from = current.cluster.ids, to = new.cluster.ids)

metadata_filt_DC_N_T <- metadata_filt_DC_N_T[metadata_filt_DC_N_T$cancer %in% c("BRCA","CRC","ESCA","STAD","HCC","LUNG","OV-FTC","PAAD","KIDNEY","THCA","UCEC"),]
cell_num <- ddply(metadata_filt_DC_N_T, .(patient,cancer,tissue), nrow)
colnames(cell_num) <- c("patient","cancer","tissue","Total")

metadata_filt_DC_N_T_LAMP3 <- metadata_filt_DC_N_T[metadata_filt_DC_N_T$DC_type %in% c("cDC3"),]
cell_num_cluster <- ddply(metadata_filt_DC_N_T_LAMP3, .(patient,cancer,tissue), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","tissue","number")
cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer~tissue,value.var=c('Proportion'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer"))

p <- ggpaired(cell_num_cluster_summary_formated_melt, x = "variable", y = "value",
          color='variable',palette = "jco", 
          line.color = "gray", line.size = 0.4,
           short.panel.labs = TRUE) +facet_wrap(~cancer, scales = "free",strip.position="top",ncol =3)+ ylab("Proportion of LAMP3+ cDCs in cDCs")+ xlab("")

p + stat_compare_means(paired = TRUE,method.args = list(alternative = "less"))


### plot macrophage subset distribution for lung cancer and UCEC
## macrophage from Lung cancer and UCEC
metadata_filt_N_T_cancer <- metadata_filt_N_T[metadata_filt_N_T$cancer %in% c("Lung","UCEC","THCA"),]
metadata_filt_N_T_cancer$tissue <- factor(metadata_filt_N_T_cancer$tissue)
cell_num <- ddply(metadata_filt_N_T_cancer, .(patient,cancer,tissue), nrow)
colnames(cell_num) <- c("patient","cancer","tissue","Total")

metadata_filt_SPP1 <- metadata_filt_N_T_cancer[metadata_filt_N_T_cancer$MajorCluster %in% c("M10_Macro_SPP1","M09_Macro_SPP1"),]
cell_num_cluster <- ddply(metadata_filt_SPP1, .(patient,cancer,tissue), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","tissue","number")
cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary$Proportion <- round((cell_num_cluster_summary$number/cell_num_cluster_summary$Total)*100, 2)

## paired wilcox test
library(ggpubr)
cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer~tissue,value.var=c('Proportion'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer"))

p <- ggpaired(cell_num_cluster_summary_formated_melt, x = "variable", y = "value",
          color='variable',palette = "jco", 
          line.color = "gray", line.size = 0.4,
           short.panel.labs = TRUE) +facet_wrap(~cancer, scales = "free",strip.position="top",ncol =3)+ ylab("Proportion of Macro_SPP1")+ xlab("")
p + stat_compare_means(paired = TRUE)



########### only consider composition of myeloid cell in tumor tissues from different tumor types
setwd("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate")
metadata <- read.csv("metadata.csv")
metadata_filt <- metadata[metadata$Cluster!='Myeloid',]
metadata_filt <- metadata_filt[metadata_filt$tissue=='T',]

current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
metadata_filt$cancer <- plyr::mapvalues(x = metadata_filt$cancer, from = current.cluster.ids, to = new.cluster.ids)


cell_num <- ddply(metadata_filt, .(patient,cancer), nrow)
colnames(cell_num) <- c("patient","cancer","Total")

cell_num_cluster <- ddply(metadata_filt, .(patient,cancer,Cluster), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","cluster","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+cancer+Total~cluster,value.var=c('number'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","cancer","Total"))

cell_num_cluster_summary_formated_melt <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$Total >=100,]
cell_num_cluster_summary_formated_melt$Proportion <- round((cell_num_cluster_summary_formated_melt$value/cell_num_cluster_summary_formated_melt$Total)*100, 2)

# Color panel -----------------
  c54 <- c('dodgerblue2','green4','#E31A1C','#6A3D9A','#FF7F00',
         '#FB9A99','#CAB2D6','khaki2','deeppink1','blue1',      
         'steelblue4','green1','yellow4','yellow3','forestgreen',
         'red2','orange','cornflowerblue', 'magenta','darkolivegreen4',
         'indianred1','tan4','darkblue','mediumorchid1','firebrick4',
         'yellowgreen','lightsalmon','tan3','tan1','darkgray',
         'wheat4','#DDAD4B','chartreuse','seagreen1','moccasin',
         'mediumvioletred','seagreen','cadetblue1','darkolivegreen1','tan2',
         'tomato3','#7CE3D8','gainsboro','gold1','skyblue2',
         'palegreen2','#FDBF6F','gray70','darkorange4','orchid1',
         'darkturquoise','maroon','brown','black')


c54 <- c("BRCA" = 'dodgerblue2',"ESCA"='green4',"STAD"='#E31A1C',"LUNG"='#6A3D9A',"OV-FTC"='#FF7F00',
         "CRC"='#FB9A99',"PAAD"='#CAB2D6',"KIDNEY"='khaki2',"THCA"='deeppink1',"UCEC"='blue1',      
         "HCC"='steelblue4',"NPC"='green1',"MEL"='yellow4',"MYE"='yellow3',"LYM"='forestgreen')

type <- "Mast"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p1 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 70, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
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

type <- "pDC"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p2 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 50, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
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
type <- "cDC"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p3 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 50, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
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
	  
	  
type <- "Mo/Mq"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p4 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
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

library(scater)
pdf("/data2/csj/Pan_Myeloid/A20191105/P1_fraction.pdf",width=4.78,height=6.23)
multiplot(p1,p2,cols=1)
dev.off()
pdf("/data2/csj/Pan_Myeloid/A20191105/P2_fraction.pdf",width=4.78,height=6.23)
multiplot(p3,p4,cols=1)
dev.off()

multiplot(p4,p2,p3,p1,cols=1)
	  

#### only consider DC cells
## DC subset comparison
## DC subset cluster distribution
metadata_filt_DC <- metadata_filt[metadata_filt$Cluster %in% c("cDC"),]
metadata_filt_DC$DC_type <- unlist(lapply(strsplit(metadata_filt_DC$MajorCluster,"_"),function(x){x[2]}))

cell_num <- ddply(metadata_filt_DC, .(patient,cancer), nrow)
colnames(cell_num) <- c("patient","cancer","Total")

cell_num_cluster <- ddply(metadata_filt_DC, .(patient,cancer,DC_type), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","cluster","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0
cell_num_cluster_summary <- cell_num_cluster_summary[cell_num_cluster_summary$Total >=50,]

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+Total+cancer~cluster,value.var=c('number'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","Total","cancer"))

cell_num_cluster_summary_formated_melt$Proportion <- round((cell_num_cluster_summary_formated_melt$value/cell_num_cluster_summary_formated_melt$Total)*100, 2)

type <- "cDC1"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p1 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
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
type <- "cDC2"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p2 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
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

type <- "cDC3"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

p3 <- ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of ", type," in tumor"))+ theme(
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
library(scater)
multiplot(p1,p2,p3,cols=1)
pdf("/data2/csj/Pan_Myeloid/A20191105/cDC_fraction.pdf",width=4.78,height=9)
multiplot(p1,p2,p3,cols=1)
dev.off()
  
########### only consider composition of myeloid cell in tumor tissues from different tumor types and calcute anogigenesis and phagocytosis ratio

setwd("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate")
metadata <- read.csv("metadata.csv")
metadata_filt <- metadata[metadata$Cluster!='Myeloid',]
metadata_filt <- metadata_filt[metadata_filt$tissue=='T',]

current.cluster.ids <- c("BC","ESCA","Gastric","Lung","OV","CRC","PACA","RC","THCA","UCEC","HCC","NPC","Mela","MM","L")
new.cluster.ids <- c("BRCA","ESCA","STAD","LUNG","OV-FTC","CRC","PAAD","KIDNEY","THCA","UCEC","HCC","NPC","MEL","MYE","LYM")
metadata_filt$cancer <- plyr::mapvalues(x = metadata_filt$cancer, from = current.cluster.ids, to = new.cluster.ids)

metadata_filt$class <- "other"
metadata_filt[metadata_filt$cancer=='BRCA' & metadata_filt$MajorCluster =='M09_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='BRCA' & metadata_filt$MajorCluster =='M08_Macro_CLEC5A',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='BRCA' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='OV-FTC' & metadata_filt$MajorCluster =='M08_Macro_MARCO',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='OV-FTC' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='PAAD' & metadata_filt$MajorCluster =='M08_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='PAAD' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='KIDNEY' & metadata_filt$MajorCluster =='M08_Macro_FN1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='KIDNEY' & metadata_filt$MajorCluster =='M11_Macro_SLC40A1',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='STAD' & metadata_filt$MajorCluster =='M08_Macro_IL1B',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='STAD' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='THCA' & metadata_filt$MajorCluster =='M09_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='THCA' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='UCEC' & metadata_filt$MajorCluster =='M10_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='UCEC' & metadata_filt$MajorCluster =='M11_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='CRC' & metadata_filt$MajorCluster =='M11_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='CRC' & metadata_filt$MajorCluster =='M12_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='ESCA' & metadata_filt$MajorCluster =='M09_Macro_IDO1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='ESCA' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='HCC' & metadata_filt$MajorCluster =='M08_Macro_FCN1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='HCC' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='LUNG' & metadata_filt$MajorCluster =='M09_Macro_SPP1',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='LUNG' & metadata_filt$MajorCluster =='M10_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='LYM' & metadata_filt$MajorCluster =='M06_Macro_IL1B',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='LYM' & metadata_filt$MajorCluster =='M07_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='MEL' & metadata_filt$MajorCluster =='M07_Macro_VCAN',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='MEL' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'

metadata_filt[metadata_filt$cancer=='MYE' & metadata_filt$MajorCluster =='M08_Macro_IL1B',]$class = 'anogigenesis'
metadata_filt[metadata_filt$cancer=='MYE' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'phagocytosis'


metadata_filt[metadata_filt$cancer=='NPC' & metadata_filt$MajorCluster =='M09_Macro_C1QA',]$class = 'anogigenesis'

cell_num_cluster <- ddply(metadata_filt, .(patient,cancer,class), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","cluster","number")

cell_num_cluster_formated <- dcast(cell_num_cluster,patient+cancer~cluster,value.var=c('number'))
cell_num_cluster_formated[is.na(cell_num_cluster_formated)] <- 0
cell_num_cluster_formated$Total <- cell_num_cluster_formated$anogigenesis + cell_num_cluster_formated$phagocytosis + cell_num_cluster_formated$other
cell_num_cluster_formated <- cell_num_cluster_formated[cell_num_cluster_formated$Total >=100,]
cell_num_cluster_formated$Ratio <- round(log2(cell_num_cluster_formated$anogigenesis/cell_num_cluster_formated$phagocytosis),2)


median_table <- ddply(cell_num_cluster_formated,.(cancer), function(x){median(x$Ratio)})
median_table_sorted <- median_table[order(median_table$V1),]
cell_num_cluster_formated$cancer <- factor(cell_num_cluster_formated$cancer, levels=median_table_sorted$cancer)


ggboxplot(cell_num_cluster_formated, x = "cancer", y = "Ratio",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 6, label.x = 2) + ggtitle(paste0("Anogigenesis/phagocytosis macrophages"))+ theme(
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
      )
	  
## calcutate the anogigenesis fraction in all myeloid cell
cell_num <- ddply(metadata_filt, .(patient,cancer), nrow)
colnames(cell_num) <- c("patient","cancer","Total")

cell_num_cluster <- ddply(metadata_filt, .(patient,cancer,class), nrow)
colnames(cell_num_cluster) <- c("patient","cancer","cluster","number")

cell_num_cluster_summary <-merge(cell_num,cell_num_cluster, all.x = TRUE)
cell_num_cluster_summary[is.na(cell_num_cluster_summary)] <- 0

cell_num_cluster_summary_formated <- dcast(cell_num_cluster_summary,patient+Total+cancer~cluster,value.var=c('number'))
cell_num_cluster_summary_formated[is.na(cell_num_cluster_summary_formated)] <- 0
cell_num_cluster_summary_formated_melt <- melt(cell_num_cluster_summary_formated,id=c("patient","Total","cancer"))

cell_num_cluster_summary_formated_melt <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$Total >=100,]
cell_num_cluster_summary_formated_melt$Proportion <- round((cell_num_cluster_summary_formated_melt$value/cell_num_cluster_summary_formated_melt$Total)*100, 2)


type <- "anogigenesis"
df <- cell_num_cluster_summary_formated_melt[cell_num_cluster_summary_formated_melt$variable == type,]
median_table <- ddply(df,.(cancer), function(x){median(x$Proportion)})
median_table_sorted <- median_table[order(median_table$V1),]
df$cancer <- factor(df$cancer, levels=median_table_sorted$cancer)

pdf("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate/angiogenesis-associated.pdf",width=6.19,height=3.51)
ggboxplot(df, x = "cancer", y = "Proportion",
                color = "cancer", palette =c54,
                add = "jitter") + stat_compare_means(label.y = 100, label.x = 2) + ggtitle(paste0("Proportion of angiogenesis-associated macrophages"))+ theme(
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
      )
dev.off()	  