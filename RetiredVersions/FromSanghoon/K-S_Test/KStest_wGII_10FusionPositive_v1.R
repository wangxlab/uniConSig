library(data.table)
library(dplyr)
#################################################################################
##required for all calculations
setwd("/zfs2/xiaosongwang/sal170/4-1_BCL2L14-ETV6_Cibersort_wGII_KStest/5_wGII_fusionPositive_KStest")
source("CSEA.modulesVb2.5.R")

gmtfile="FusionType10_TCGAIDs.gmt"
weightFile<-"TCGAID_wGII_weight.txt"
############################################################################################

feature.list=read_gmt(gmtfile,min=0)
weightData<-fread(weightFile,stringsAsFactors=F,sep="\t")
weight=weightData$wGII
names(weight)=weightData$TCGAID

ks.result=run.weightedKS(weight,signed=F,feature.list,minsize=0,correct.overfit=FALSE,correct.outlier=FALSE,transformNegWeight=FALSE)
	