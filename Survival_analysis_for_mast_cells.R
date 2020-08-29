## Mast cell survival annalysis
libary(survival)
library(survminer)

genelist <- c("KIT","TPSAB1","CPA3","TPSB2")

cancers <- c("LUNG","KIDNEY","BRCA","STAD","SKCM","LIHC","CRC","ESCA","PAAD","THCA")
res_list <- list()
for(cancerType in cancers){
#cancerType <- "COAD"
res <- DoSurvivalPlot(cancer = cancerType, genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("gender", "stage", "age"), survival.used = "OS", signature.name = "GeneSet")
res_list[[cancerType]] <- res
}

cancers <- c("OV","UCEC")
for(cancerType in cancers){
res <- DoSurvivalPlot(cancer = cancerType, genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("age"), survival.used = "OS", signature.name = "GeneSet")
res_list[[cancerType]] <- res
}

### summary all result
cancers <- c("LUNG","KIDNEY","BRCA","STAD","SKCM","LIHC","CRC","ESCA","PAAD","THCA","OV","UCEC")
for(cancerType in cancers){
	p_value <- round(summary(res_list[[cancerType]]$cox_result)$coefficients[1,5],2)
	HR <- round(summary(res_list[[cancerType]]$cox_result)$coefficients[1,2],2)
	cat(cancerType," P-value:", p_value,"; HR: ", HR,"\n")
}


## PFI index

cancers <- c("LUNG","KIDNEY","BRCA","STAD","SKCM","LIHC","CRC","ESCA","PAAD","THCA")
res_PFI_list <- list()
for(cancerType in cancers){
#cancerType <- "COAD"
res <- DoSurvivalPlot(cancer = cancerType, genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("gender", "stage", "age"), survival.used = "PFI", signature.name = "GeneSet")
res_PFI_list[[cancerType]] <- res
}

cancers <- c("OV","UCEC")
for(cancerType in cancers){
res <- DoSurvivalPlot(cancer = cancerType, genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("age"), survival.used = "PFI", signature.name = "GeneSet")
res_PFI_list[[cancerType]] <- res
}

## DFI index
## SKCM without DFI
cancers <- c("LUNG","KIDNEY","BRCA","STAD","LIHC","CRC","ESCA","PAAD","THCA")
res_DFI_list <- list()
for(cancerType in cancers){
#cancerType <- "COAD"
res <- DoSurvivalPlot(cancer = cancerType, genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("gender", "stage", "age"), survival.used = "DFI", signature.name = "GeneSet")
res_DFI_list[[cancerType]] <- res
}

cancers <- c("OV","UCEC")
for(cancerType in cancers){
res <- DoSurvivalPlot(cancer = cancerType, genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("age"), survival.used = "DFI", signature.name = "GeneSet")
res_DFI_list[[cancerType]] <- res
}


## DSS index
cancers <- c("LUNG","KIDNEY","BRCA","STAD","SKCM","LIHC","CRC","ESCA","PAAD","THCA")
res_DSS_list <- list()
for(cancerType in cancers){
#cancerType <- "COAD"
res <- DoSurvivalPlot(cancer = cancerType, genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("gender", "stage", "age"), survival.used = "DSS", signature.name = "GeneSet")
res_DSS_list[[cancerType]] <- res
}

cancers <- c("OV","UCEC")
for(cancerType in cancers){
res <- DoSurvivalPlot(cancer = cancerType, genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("age"), survival.used = "DSS", signature.name = "GeneSet")
res_DSS_list[[cancerType]] <- res
}



df_HR <- data.frame()
cancers <- c("LUNG","KIDNEY","BRCA","STAD","SKCM","LIHC","CRC","ESCA","PAAD","THCA","OV","UCEC")
for(cancerType in cancers){
	if(cancerType == 'SKCM'){
		DFI_p_value <- NA
		DFI_HR <- NA
	}else{
		DFI_p_value <- round(summary(res_DFI_list[[cancerType]]$cox_result)$coefficients[1,5],2)
		DFI_HR <- round(summary(res_DFI_list[[cancerType]]$cox_result)$coefficients[1,2],2)
	}
	PFI_p_value <- round(summary(res_PFI_list[[cancerType]]$cox_result)$coefficients[1,5],2)
	PFI_HR <- round(summary(res_PFI_list[[cancerType]]$cox_result)$coefficients[1,2],2)
	OS_p_value <- round(summary(res_list[[cancerType]]$cox_result)$coefficients[1,5],2)
	OS_HR <- round(summary(res_list[[cancerType]]$cox_result)$coefficients[1,2],2)
	DSS_p_value <- round(summary(res_DSS_list[[cancerType]]$cox_result)$coefficients[1,5],2)
	DSS_HR <- round(summary(res_DSS_list[[cancerType]]$cox_result)$coefficients[1,2],2)
	df_HR <- rbind(df_HR, c(cancerType,PFI_p_value,PFI_HR,OS_p_value,OS_HR,DFI_p_value,DFI_HR,DSS_p_value,DSS_HR))
}
colnames(df_HR) <- c("cancer","PFI_Pvalue","PFI_HR","OS_Pvalue","OS_HR", "DFI_Pvalue","DFI_HR", "DSS_Pvalue","DSS_HR")

### heatmap
df_HR[10,9] <- 10^0.5
cancer_list <- df_HR$cancer
signatures <- c("PFI","OS","DFI","DSS")
m_HR <- as.matrix(t(data.frame(PFI=log10(as.numeric(df_HR$PFI_HR)),OS=log10(as.numeric(df_HR$OS_HR)),DFI=log10(as.numeric(df_HR$DFI_HR)) ,DSS=log10(as.numeric(df_HR$DSS_HR)))))
m_p <- as.matrix(t(data.frame(PFI=df_HR$PFI_Pvalue,OS=df_HR$OS_Pvalue,DFI=df_HR$DFI_Pvalue, DSS=df_HR$DSS_Pvalue)))
# colnames(m_HR) <- df_HR$cancer
# colnames(m_p) <- df_HR$cancer
conf_signif <- 0.05


suppressPackageStartupMessages(library(ggplot2))
plt=ggplot()+theme_minimal()+
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, size = 14, hjust = 1),
      axis.text.y = element_text(size = 14),
      legend.position = "right",
      plot.title = element_text(vjust = 0.5, size = 16, hjust = 0.5),
      panel.grid = element_blank()
      )+coord_fixed()+
    scale_x_continuous('',breaks=1:length(cancer_list),labels=cancer_list,expand=rep(0.1/length(cancer_list),2))+
    scale_y_continuous('',breaks=1:length(signatures),labels=signatures,expand=rep(0.1/length(signatures),2))

  plt=plt+scale_fill_gradient2(
    low = "#0080FF", high = "#FF6666", mid = "#FFFFFF",
    midpoint = 0, name="log10(HR)",limit=c(-0.6,0.6)
  )


  for(i in 1:length(signatures))plt=plt+geom_tile(aes_(1:length(cancer_list),i,fill=m_HR[i,]))
  for(i in 1:length(signatures)){
    mark=as.logical(m_p[i,]<= conf_signif)
    clr=ifelse(m_HR[i,]>0,'#ee0000','#0000ee')
    #clr=ifelse(m_HR[i,]>1,'#ee0000','#0000ee')
    if(sum(mark,na.rm = TRUE)==0)next
    plt=plt+annotate(
      'tile',x=c(1:length(cancer_list))[mark],y=i,color=clr[mark],
      fill='transparent',size=1,width=1
    )
  }
plt




DoSurvivalPlot(cancer = "UCEC", genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("age"), survival.used = "DFI", signature.name = "Mast cell signature")

DoSurvivalPlot(cancer = "LIHC", genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("gender", "stage", "age"), survival.used = "DFI", signature.name = "Mast cell signature")

DoSurvivalPlot(cancer = "PAAD", genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c("gender", "stage", "age"), survival.used = "DFI", signature.name = "Mast cell signature")

DoSurvivalPlot(cancer = "OV", genelist, normalize.by.gene = "PTPRC", low.quantile = 0.4, high.quantile = 0.6, column.used = c( "age"), survival.used = "DFI", signature.name = "Mast cell signature")

