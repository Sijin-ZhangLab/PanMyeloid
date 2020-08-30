## /share/tools/R/R-3.6.0/bin/R
.libPaths("/data2/csj/tools/Rlib3.6")

suppressPackageStartupMessages({
library("reshape2")
library("plyr")
library("dplyr")
library("ggpubr")
library("data.table")
library("MASS")
library("R.utils")
library("nlme")
#library("MuMIn")
library("forestplot")
library("glmmLasso")
#library("lmmen")
})

options(stringsAsFactors=F)
## read TMB
tmb = read.table("/data2/csj/Pan_Myeloid/A20191105/mutation/TMB.merge.txt", header=F, check.names=F, stringsAsFactors=F)
colnames(tmb) = c("patient","cancerType","TMB")
tmb$patient = gsub("-",".",tmb$patient)
tmb$patient = gsub("BC","BRCA",tmb$patient)
tmb$patient = gsub("PACA","PAAD",tmb$patient)
old.name <- c("BC", "CHOL", "CRC", "ESCA", "HCC", "LUNG", "FTC", "OV", "PACA", "RC", "THCA", "UCEC")
new.name <- c("BRCA", "CHOL", "CRC", "ESCA", "HCC", "LUNG", "FTC", "OV", "PAAD", "RC", "THCA", "UCEC")
tmb$cancerType <- plyr::mapvalues(x = tmb$cancerType, from = old.name, to = new.name)


## read cellInfo
cellInfo = read.csv("/data2/csj/Pan_Myeloid/A20191105/scanorama_integrate/metadata.csv")
cellInfo = cellInfo[cellInfo$Cluster != 'Myeloid',]

cellInfo = cellInfo[cellInfo$tissue=="T",]
cellInfo = cellInfo[cellInfo$source=="ZhangLab",] 
current.cluster.ids <- c("ESCA", "L", "MM", "OV", "PACA", "RC", "THCA", "UCEC")
new.cluster.ids <- c("ESCA", "LYM", "MM", "OV", "PAAD", "RC", "THCA", "UCEC")
cellInfo$cancer <- plyr::mapvalues(x = cellInfo$cancer, from = current.cluster.ids, to = new.cluster.ids)
cellInfo <- cellInfo[c('cancer','patient','Cluster','MajorCluster','library_id')]
cellInfo$patient <- paste0(cellInfo$cancer,".",cellInfo$patient)

cellInfo$subcluster <- unlist(lapply(strsplit(cellInfo$MajorCluster,"_"),function(x){x[2]}))
old.cluster.name <- c("Mono", "Mast", "Macro", "pDC", "cDC3", "cDC2", "cDC1", "Monolike")
new.cluster.name <- c("Mo/Mq", "Mast", "Mo/Mq", "pDC", "cDC3", "cDC2", "cDC1", "Mo/Mq")
cellInfo$subcluster <- plyr::mapvalues(x = cellInfo$subcluster, from = old.cluster.name, to = new.cluster.name)
colnames(cellInfo) <- c('cancerType','patient','Majorcluster','Cluster','library_id','meta.cluster')

# subset, calculate and then merge tables
cells = cellInfo[,c("patient","meta.cluster")] %>% table %>% as.data.frame()
cells = cells %>% group_by(patient) %>% dplyr::mutate(Percent=100*Freq/sum(Freq))
cells$patient = as.character(cells$patient)
cells$meta.cluster = as.character(cells$meta.cluster)


overlap = intersect(unique(tmb$patient), unique(cells$patient))
length(overlap)  ## 35
tmb = tmb[tmb$patient %in% overlap, ]
cellInfo = cellInfo[cellInfo$patient %in% overlap, ]
cells = cells[cells$patient %in% overlap,]

dat = merge(cells, tmb[,c("patient","cancerType","TMB")], by="patient")
# TMB
## overall  
pdf("/data2/csj/Pan_Myeloid/A20191105/01_revised_data/Overall_TMB_cell_proportion.pdf",width=6.4,height=4.57)
ggplot(dat, aes(x=TMB, y=Percent)) + 
  geom_point() +
  stat_cor(size=3,method = "spearman") +
  geom_smooth(method='lm',formula=y~x) +
  theme_bw() + 
  facet_wrap(.~meta.cluster, scales="free", nrow=2)+
  theme_bw()+ theme(strip.text.x = element_text(size = 12, colour = "red"),
		panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

dev.off()

# Spec. Mut.
## load files, make the dataframe
mut = read.table("/data2/csj/Pan_Myeloid/A20191105/mutation/mut_mtx.flt.fmt.txt", sep="\t", header=T, check.names=F, stringsAsFactors=F)
mut$Gene = gsub("-","_", mut$Gene)
colnames(mut) = gsub("-",".",colnames(mut))
colnames(mut) = gsub("PACA","PAAD",colnames(mut))
colnames(mut) = gsub("BC","BRCA",colnames(mut))
rownames(mut) = mut$Gene
#
flag = mut[,overlap]
flag = ifelse(flag==".", 0, 1)
mut$Num2 = rowSums(flag)
mut$Num2r = round(mut$Num2/length(overlap), 2)
#
census = fread("/data2/csj/Pan_Myeloid/A20191105/mutation/Census_all_Wed_Oct_30_07_30_28_2019.GRCh37_COSMICv90.tsv", sep="\t", header=T, stringsAsFactors=F)
census$`Gene Symbol` = gsub("-","_", census$`Gene Symbol`)

## choose gene used
# set1: Num2r>=0.1 && in Census table
table(mut$Num2)
genes = mut[mut$Num2r >= 0.1,"Gene"]
genes = genes[genes %in% census$`Gene Symbol`]
genes = genes[!grepl("[;,]",genes)]
length(genes)
print(genes)

#reshape
mut = mut[genes, overlap]
mut = as.data.frame(t(ifelse(mut==".", 0, 1)))

dat2 = dcast(cells, patient ~ meta.cluster, value.var="Percent")
dat2[is.na(dat2)] = 0
dat2 = merge(unique(dat[,c("patient","cancerType")]), dat2, by="patient")
clusters = colnames(dat2)[3:ncol(dat2)]

mut$patient = rownames(mut)
dat2 = merge(dat2, mut, by="patient")
dat2 = merge(dat2, tmb[,c("patient","TMB")], by="patient")
dat2$cancerType = factor(dat2$cancerType)

## used save lm result
allRes <<- as.data.frame(matrix(NA, nrow=0, ncol=5))
colnames(allRes) = c("Gene", "Cluster", "Estimate", "EstimateCI", "Pvalue")

## define functions
## lasso variables selection, lambda choosed using 3-fold cv
runLasso = function(par.dat2, par.tarClu, par.genes){
  # https://stats.stackexchange.com/questions/70753/correct-estimation-of-arguments-for-glmmlasso-function/134559
  set.seed(123)
  N = nrow(par.dat2)
  ind = sample(N,N)
  lambda = seq(500, 0, by=-5)
  family = gaussian(link="identity")
  #
  par.dat2$y = par.dat2[,par.tarClu[1]]
  if (length(par.tarClu) > 1){
    for(i in par.tarClu[2:length(par.tarClu)]){
       par.dat2$y = par.dat2$y + par.dat2[,i]
    }
  }
  formula = as.formula(paste("y ~", paste(c("TMB",par.genes), collapse='+')))
  #
  kk = 3
  nk = floor(N/kk)
  Devianz_ma = matrix(Inf,ncol=kk,nrow=length(lambda))
  ## first fit good starting model
  PQL = glmmPQL(y~1, random=~1|cancerType, family=family, data=par.dat2, verbose=F)
  Delta.start = c(as.numeric(PQL$coef$fixed), rep(0,length(par.genes)+1), as.numeric(t(PQL$coef$random$cancerType)))
  Q.start = as.numeric(VarCorr(PQL)[1,1])
  ## train and test
  for(j in 1:length(lambda)){
    for (i in 1:kk){
      if (i < kk){
        indi = ind[(i-1)*nk+(1:nk)]
      }else{
        indi = ind[((i-1)*nk+1):N]
      }
      
      par.dat2.train = par.dat2[-indi,]
      par.dat2.test = par.dat2[indi,]
  
      glm2 = try(glmmLasso(formula, rnd=list(cancerType=~1), family=family, data=par.dat2.train, lambda=lambda[j], 
                  switch.NR=F, final.re=TRUE, control=list(start=Delta.start, q_start=Q.start)), silent=TRUE) 
        
      if(class(glm2)!="try-error"){  
        y.hat = predict(glm2, par.dat2.test)    
        Devianz_ma[j,i] = sum(family$dev.resids(par.dat2.test$y, y.hat, wt=rep(1,length(y.hat))))
      }
    }
    #print(paste("Iteration ", j, sep=""))
    #print(sum(Devianz_ma[j,]))
  }
  ## get best lambda
  Devianz_vec = apply(Devianz_ma,1,sum)
  opt2 = which.min(Devianz_vec)
  ## final lasso regression
  glm2_final = glmmLasso(formula, rnd=list(cancerType=~1), family=family, data=par.dat2, lambda=lambda[opt2],
                         switch.NR=F, final.re=TRUE, control=list(start=Delta.start,q_start=Q.start))
  
  return(glm2_final)
}

## show LASSO model info
showLasso = function(par.dat2, par.tarClu, par.lm){
  sum = summary(par.lm)
  
  tab = as.data.frame(sum$coefficients)
  tab = tab[tab$Estimate!=0,]
  tab = tab[!grepl("^(TMB|\\(Intercept\\))$",rownames(tab)),]
  #
  if (nrow(tab)==0){
    print ("No good model")
    return(NULL)
  }
  print(sum)
  #
  tab$Gene = rownames(tab)
  tab = tab[order(tab$Estimate, decreasing=T), ]
  #
  tab$upper = tab$Estimate + 2*tab$StdErr
  tab$lower = tab$Estimate - 2*tab$StdErr
  tab$sig = ifelse(tab$p.value <0.05, "*", " ")
  #tab$sig = ifelse(tab$p.value <0.01, "**", tab$sig)
  tab$Pval = sprintf("%s%s", tab$p.value, tab$sig)
  tab$Info = sprintf("%.2f (%.2f - %.2f)%s", tab$Estimate, tab$lower, tab$upper, tab$sig)
  tab$Estimate = sprintf("%.2f", tab$Estimate)
  tab$StdErr = sprintf("%.2f", tab$StdErr)

  # forestplot
  print(
    forestplot(rbind(c(c("Gene", NA, "Esitimate (95% CI)")), cbind(tab$Gene, NA, tab$Info)),
      c(NA,tab$Estimate), c(NA,tab$lower), c(NA,tab$upper),
      graph.pos=2, xlab="Esitimate", is.summary=c(T, rep(F,nrow(tab))), boxsize=0.15, 
      hrzl_lines=T, new_page=T,
      col=fpColors(box=c("#1c61b6","black"), lines="#1c61b6", zero="gray50"), 
      lwd.ci=2, ci.vertices=T, ci.vertices.height = 0.05,
      txt_gp=fpTxtGp(label=gpar(cex=1), ticks=gpar(cex=0.7), xlab=gpar(cex=1), title=gpar(cex=1)),
      xticks=seq(from=floor(min(c(0,tab$lower))), length.out=10,
                  by=round((ceiling(max(c(0,tab$upper)))-floor(min(c(0,tab$lower))))/10, 2))
    )
  )
  # save and return
  tab2save(par.dat2, par.tarClu, tab, allRes)
  return(tab$Gene)
}

## boxplot for selected genes
loop.box = function(par.dat2, par.tarClu, par.genes){
  for (i in par.genes){
    p = ggboxplot(par.dat2, x=i, y=paste(tarClu,collapse="+"), color=i, add="jitter", palette="jco", outlier.shape=NA) + 
      stat_compare_means(comparisons=list(c("0","1"))) + 
      theme(legend.position="right", plot.title=element_text(size=9, vjust = -1)) +
      ylab("Percent") + xlab("Mutation") 
    print(p)
  }
}

## save res
tab2save = function(par2.dat2, par2.tarClu, par2.tab, par2.allRes){
  # cid
  cid = par2.tarClu
  if (length(par2.tarClu) > 1){
    cid = paste(par2.tarClu, collapse="+")
  }
  # test.data
  test.dat = rowSums(par2.dat2[,par2.tarClu,drop=F])
  # calculate
  tmp = par2.tab
  tmp$Cluster = cid
  tmp$EstimateCI = gsub("\\*$","",tmp$Info)
  tmp$EstimateCI = sub(" - ", " ~ ",tmp$EstimateCI)
  tmp$Pvalue = NA
  for (i in tmp$Gene){
    t = wilcox.test(test.dat[ par2.dat2[,i]==0 ], test.dat[ par2.dat2[,i]==1])
    tmp[tmp$Gene==i,"Pvalue"] = sprintf("%.3f",t$p.value)
  }
  # save
  this.res = tmp[,c("Gene", "Cluster", "Estimate", "EstimateCI", "Pvalue")]
  allRes <<- unique(rbind(par2.allRes, this.res))
}



## Lasso
tarClu = c("cDC1")
lm.lasso = runLasso(dat2, tarClu, genes)
lm.genes = showLasso(dat2, tarClu, lm.lasso)
print(lm.genes)
#loop.box(dat2, "cDC3", "MUC4")


tarClu = c("cDC2")
lm.lasso = runLasso(dat2, tarClu, genes)
lm.genes = showLasso(dat2, tarClu, lm.lasso)
print(lm.genes)

tarClu = c("cDC3")
lm.lasso = runLasso(dat2, tarClu, genes)
lm.genes = showLasso(dat2, tarClu, lm.lasso)
pdf("cDC3_MUC4.pdf",width=3,height=3)
loop.box(dat2, tarClu, "MUC4")
dev.off()

### plot with cancer type
# ggboxplot(dat2, x='MUC4', y=paste('cDC3',collapse="+"), add="jitter", palette="jco", outlier.shape=NA,add.params = list(color = 'cancerType'))+ 
      # stat_compare_means(comparisons=list(c("0","1"))) + 
      # theme(legend.position="right", plot.title=element_text(size=9, vjust = -1)) +
      # ylab("The proportion of LAMP3+ DC") + xlab("MUC4 Mutation")

tarClu = c("pDC")
lm.lasso = runLasso(dat2, tarClu, genes)
lm.genes = showLasso(dat2, tarClu, lm.lasso)
print(lm.genes)

tarClu = c("Mast")
lm.lasso = runLasso(dat2, tarClu, genes)
lm.genes = showLasso(dat2, tarClu, lm.lasso)
print(lm.genes)

tarClu = c("Mo/Mq")
lm.lasso = runLasso(dat2, tarClu, genes)
lm.genes = showLasso(dat2, tarClu, lm.lasso)
print(lm.genes)


## Final save result

allRes$Estimate.p = as.numeric(allRes$Estimate)
allRes$Estimate.p = ifelse(allRes$Estimate.p > 5, 5, allRes$Estimate.p)
allRes$Estimate.p = ifelse(allRes$Estimate.p < -5, -5, allRes$Estimate.p)
#
allRes$EstimateCI.p = sub(" ","\n", allRes$EstimateCI)
allRes$EstimateCI.p = gsub(" ","", allRes$EstimateCI.p)
allRes$EstimateCI.p = gsub("(\\d\\.\\d)\\d","\\1", allRes$EstimateCI.p)
#
allRes$sig = ifelse(allRes$Pvalue<0.05, "P<0.05", NA)
allRes$sig = factor(allRes$sig, levels=c("P<0.05"))
sig.color = c("#EE2C2C")
names(sig.color) = c("P<0.05")
#
pdf("mutation_myeloid_subset_new.pdf", width=5, height=2)
ggplot(allRes, aes(x=Gene, y=Cluster, fill=Estimate.p)) +
                geom_tile(size=3) +
                geom_tile(aes(color=sig),size=1.2, alpha=0) +
                geom_text(aes(label=Estimate), angle=0, color="black", size=3) +
                scale_fill_gradient2( high="#EEAEEE", mid="#FFFFFF", low="#1C86EE", midpoint=0, limits=c(-5,5), name="Coef") +
                scale_colour_manual(name="significance",values=sig.color) +
                theme_bw() + 
                theme(axis.text.x=element_text(size=8, vjust=1, hjust=0, angle=300)) +
                xlab('') + ylab("cellType")+theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))

dev.off()				



########################################################################################################3
## pathway level regression
# Spec. Mut.
## load files, make the dataframe
mut = read.table("/data2/csj/Pan_Myeloid/A20191105/mutation/hallmark_mtx.txt", sep="\t", header=T, check.names=F, stringsAsFactors=F)
colnames(mut) = gsub("-",".",colnames(mut))
colnames(mut) = gsub("PACA","PAAD",colnames(mut))
colnames(mut) = gsub("BC","BRCA",colnames(mut))
colnames(mut) = gsub("PATH","Gene",colnames(mut))
rownames(mut) = mut$Gene
#
flag = mut[,overlap]
flag = ifelse(flag==".", 0, 1)
mut$Num2 = rowSums(flag)
mut$Num2r = round(mut$Num2/length(overlap), 2)

## choose pathway used
# set1: Num2r>=0.1 
table(mut$Num2)
pathways = mut[mut$Num2r >= 0.1,"Gene"]

length(pathways)
print(pathways)

#reshape
mut = mut[pathways, overlap]
mut = as.data.frame(t(ifelse(mut==".", 0, 1)))

dat2 = dcast(cells, patient ~ meta.cluster, value.var="Percent")
dat2[is.na(dat2)] = 0
dat2 = merge(unique(dat[,c("patient","cancerType")]), dat2, by="patient")
clusters = colnames(dat2)[3:ncol(dat2)]

mut$patient = rownames(mut)
dat2 = merge(dat2, mut, by="patient")
dat2 = merge(dat2, tmb[,c("patient","TMB")], by="patient")
dat2$cancerType = factor(dat2$cancerType)

## used save lm result
allRes <<- as.data.frame(matrix(NA, nrow=0, ncol=5))
colnames(allRes) = c("Gene", "Cluster", "Estimate", "EstimateCI", "Pvalue")	

## Lasso
tarClu = c("cDC1")
lm.lasso = runLasso(dat2, tarClu, pathways)
lm.pathways = showLasso(dat2, tarClu, lm.lasso)
print(lm.pathways)

tarClu = c("cDC2")
lm.lasso = runLasso(dat2, tarClu, pathways)
lm.pathways = showLasso(dat2, tarClu, lm.lasso)
print(lm.pathways)

tarClu = c("cDC3")
lm.lasso = runLasso(dat2, tarClu, pathways)
lm.pathways = showLasso(dat2, tarClu, lm.lasso)
print(lm.pathways)
pdf("Path1.pdf",width=4,height=3)
loop.box(dat2, tarClu, "HALLMARK_APOPTOSIS")
dev.off()
pdf("Path2.pdf",width=5,height=3)
loop.box(dat2, tarClu, "HALLMARK_INFLAMMATORY_RESPONSE")
dev.off()

### plot with cancer type
# ggboxplot(dat2, x='HALLMARK_APOPTOSIS', y=paste('cDC3',collapse="+"), add="jitter", palette="jco", outlier.shape=NA,add.params = list(color = 'cancerType'))+ 
      # stat_compare_means(comparisons=list(c("0","1"))) + 
      # theme(legend.position="right", plot.title=element_text(size=9, vjust = -1)) +
      # ylab("The proportion of LAMP3+ DC") + xlab("Apoptosis Mutation")
	  
# ggboxplot(dat2, x='HALLMARK_INFLAMMATORY_RESPONSE', y=paste('cDC3',collapse="+"), add="jitter", palette="jco", outlier.shape=NA,add.params = list(color = 'cancerType'))+ 
      # stat_compare_means(comparisons=list(c("0","1"))) + 
      # theme(legend.position="right", plot.title=element_text(size=9, vjust = -1)) +
      # ylab("The proportion of LAMP3+ DC") + xlab("Inflammatory Response Mutation")
	 	  
	  
tarClu = c("pDC")
lm.lasso = runLasso(dat2, tarClu, pathways)
lm.pathways = showLasso(dat2, tarClu, lm.lasso)
print(lm.pathways)

tarClu = c("Mast")
lm.lasso = runLasso(dat2, tarClu, pathways)
lm.pathways = showLasso(dat2, tarClu, lm.lasso)
print(lm.pathways)

tarClu = c("Mo/Mq")
lm.lasso = runLasso(dat2, tarClu, pathways)
lm.pathways = showLasso(dat2, tarClu, lm.lasso)
print(lm.pathways)


## Final save result

allRes$Estimate.p = as.numeric(allRes$Estimate)
allRes$Estimate.p = ifelse(allRes$Estimate.p > 5, 5, allRes$Estimate.p)
allRes$Estimate.p = ifelse(allRes$Estimate.p < -5, -5, allRes$Estimate.p)
#
allRes$EstimateCI.p = sub(" ","\n", allRes$EstimateCI)
allRes$EstimateCI.p = gsub(" ","", allRes$EstimateCI.p)
allRes$EstimateCI.p = gsub("(\\d\\.\\d)\\d","\\1", allRes$EstimateCI.p)
#
allRes$sig = ifelse(allRes$Pvalue<0.05, "P<0.05", NA)
allRes$sig = factor(allRes$sig, levels=c("P<0.05"))
sig.color = c("#EE2C2C")
names(sig.color) = c("P<0.05")
#
allRes$Gene = gsub("HALLMARK_","",allRes$Gene)

pdf("pathway_myeloid_subset_new.pdf",width=5,height=2)
ggplot(allRes, aes(x=Gene, y=Cluster, fill=Estimate.p)) +
                geom_tile(size=3) +
                geom_tile(aes(color=sig),size=1.2, alpha=0) +
                geom_text(aes(label=Estimate), angle=0, color="black", size=3) +
                scale_fill_gradient2( high="#EEAEEE", mid="#FFFFFF", low="#1C86EE", midpoint=0, limits=c(-5,5), name="Coef") +
                scale_colour_manual(name="significance",values=sig.color) +
                theme_bw() + 
                theme(axis.text.x=element_text(size=8, vjust=1, hjust=0, angle=340)) +
                xlab('') + ylab("cellType")+theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))
				
dev.off()
				
