# R-GSEA 1.1 -- Gene Set Enrichment Analysis / Broad Institute 
# Currentyly, R-GSEA 1.1 was tested on R version 3.6.0

library(data.table)
library(dplyr)
library("utils")
library("tools")
#####################################################################
##### Set working directory and call source files
setwd("/zfs2/xiaosongwang/sal170/18_uniConSig_CSEA2/4-1_GSEA_R_test/R-GSEA_Disambiguation_v1.2")
## Note, cls file should be made in 'space' delimited. "\t" delimited doesn't work. If you get this error: Error in C[[2]] : subscript out of bounds
## R source program 
source("CSEA.modulesV2.3.R")
source("R-GSEA.1.2_GeneRanking_SH.R", verbose=T, max.deparse.length=9999)

#####################################################################
##### These are for GSEA
input.ds<-"../../scDataset/scData_quiescentVsActive_noDupRow.gct"	# Input gene expression dataset file in GCT format
input.cls<-"../../scDataset/scData_quiescentVsActive.cls"		# Input class vector (phenotype) file in CLS format. Tab delimited works
gmtfile="../../ConceptDb/ConceptDb20190624.rmCa.gmt"
compare.list=c(read_gmt("../../PathwayDb/h.all.v6.2.symbols.gmt",min=10),read_gmt("../../PathwayDb/c2.cp.v6.2.symbols.gmt",min=10))
##################################################

## Run GSEA_SignalToNoise: it takes about 1~2 minutes
SignalToNoise<-GSEA_SignalToNoise(   ## Input/Output Files :-------------------------------------------
	input.ds=input.ds,	# Input gene expression dataset file in GCT format
	input.cls=input.cls		# Input class vector (phenotype) file in CLS format. Tab delimited doesn't work
)
#---------------------------

#########################################################
##required for all calculations
feature.list=read_gmt(gmtfile,min=10)
feature.preCalfile=paste(gsub(".gmt","",gmtfile),".Ochiai.min10.preCal.gmt",sep="")
preCalmatrix<- as.matrix(unlist(readLines(feature.preCalfile)[-1]),col=1)

ks.result=run.weightedKS(weight=SignalToNoise,signed=T,feature.list,minsize=10,correct.overfit=FALSE,correct.outlier=TRUE)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list,outfile,p.cut=0.01,q.cut=0.25,NES.cut=0,power=1,root=1,ECNpenalty=0.5,correct.overfit=FALSE)

up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)

up.disambiguate<-disambiguation(CSEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,topn=min(c(30,nrow(up.CSEA.result))),p.cut=0.01)
down.disambiguate<-disambiguation(CSEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,topn=min(c(30,nrow(down.CSEA.result))),p.cut=0.01)



