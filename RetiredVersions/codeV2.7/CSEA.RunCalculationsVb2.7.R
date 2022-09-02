##IMPORTANT: this code block is required for all calculations including both uniConSig and CSEA.
setwd("R:/CaGenome/Pipelines/DrWangPipeline/CSEA2/")
source("codeV2.7/CSEA.modulesVb2.7.R")
gmtfile="ConceptDb/ConceptDb20190624.rmCa.gmt"
feature.list=read_concepts(gmtfile,min=10)
feature.preCalfile=paste(gsub(".gmt","",gmtfile),".Ochiai.min10.preCal.gmt",sep="")
#batch_calSimilarity(feature.list=feature.list,feature.preCalfile=feature.preCalfile,cutoff=0.05,minsize=10,method="Ochiai")
preCalmatrix<- as.matrix(unlist(readLines(feature.preCalfile)[-1]),col=1)

####################################################################################################################
##Part I. Calculate uniConSig scores to compute a functional relevance of 
##        human genes underlying a target gene list such as known cancer genes
####################################################################################################################
####CalUniConSig for Cancer Gene
CGC=read.table("GeneCensus/CancerGeneCensus_allSat Jan 30 17 31 29 2021.tsv",header=TRUE,sep="\t")
ca.uniConSig.all=cal.uniConSig.cut(target.list=CGC$Gene.Symbol,feature.list,preCalmatrix,minsize=10,rm.overfit=TRUE)
write.table(ca.uniConSig.all,file="GeneCensus/Cancer.UniConSig2.1.PANCAN.tsv",quote=F,row.names=FALSE,col.names = TRUE,sep="\t")
#optional: benchmark cancer gene prioritization results
benchmark=foldCrossVal.benchmarkUniConSig(CGC$Gene.Symbol,feaure.list,preCalmatrix,fold=5,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai")

####CalUniConSig for Immune Gene
immune=read.table("GeneCensus/IMMUNE_SYSTEM_PROCESS.tsv",skip=2,sep="\t")
im.uniConSig.all=cal.uniConSig.cut(target.list=immune$V1,feature.list,preCalmatrix,minsize=10,rm.overfit=TRUE)
write.table(im.uniConSig.all,file="GeneCensus/Immune.UniConSig2.1.tsv",quote=F,row.names=FALSE,col.names = TRUE,sep="\t")

####################################################################################################################
##Part II perform CSEA for pathway enrichment analysis
####################################################################################################################
###Perform D-CSEA for dichotomous experimental gene list from scDataset
compare.list=c(read_concepts("PathwayDb/h.all.v7.5.1.symbols.gmt",min=10),read_concepts("PathwayDb/c2.cp.v7.5.1.symbols.gmt",min=10))
target.list<-read.table ("scDataset/downGenes_HSC_hgncSym.txt",header=FALSE,sep="\t")$V1
uniConSig=cal.uniConSig(target.list=target.list,feature.list=feature.list,preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai")
CSEA.result<-CSEA2(setNames(as.numeric(uniConSig$uniConSig), uniConSig$subjectID),compare.list,p.cut=0.05)
disambiguate<-disambiguation.CSEA(GEA.result=CSEA.result,uniConSig.result=uniConSig,compare.list=compare.list,topn=min(c(100,nrow(CSEA.result))),p.cut=0.01) ##upPathways use TRUE/FALSE for W-CSEA or use NULL D-CSEA
assoc <- pathwayAssociation(topPathway=disambiguate[[1]]$Compare.List[1:min(c(30,nrow(disambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
heatmap=pathway.heatmap(matrixData=assoc,clustering = TRUE,fontSize=5)

###Perform W-CSEA for pathway analysis. This method is more appliable to scRNAseq. For pathway analysis of bulk gene expression data, please use GSEA.
compare.list=c(read_concepts("PathwayDb/h.all.v7.5.1.symbols.gmt",min=10),read_concepts("PathwayDb/c2.cp.v7.5.1.symbols.gmt",min=10))
limma=limmaDGE(gctFile="./scDataset/scData_quiescentVsActive.gct",clsFile="J:/GenomeIndex/Pipeline_CSEA2/scDataset/scData_quiescentVsActive.cls")
weight=limma$Signed.Q.Value
names(weight)=row.names(limma)
ks.result=run.weightedKS(weight,signed=T,feature.list,minsize=10,correct.overfit=TRUE,correct.outlier=TRUE)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list,outfile,p.cut=0.01,q.cut=0.25,NES.cut=0,power=1,root=1,ECNpenalty=0.5,correct.overfit=FALSE)
up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)
up.disambiguate<-disambiguation.CSEA(GEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,topn=min(c(30,nrow(up.CSEA.result))),p.cut=0.01)
down.disambiguate<-disambiguation.CSEA(GEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,topn=min(c(30,nrow(down.CSEA.result))),p.cut=0.01)
up.assoc <- pathwayAssociation(topPathway=up.disambiguate[[1]]$Compare.List[1:min(c(30,nrow(up.disambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
down.assoc <- pathwayAssociation(topPathway=down.disambiguate[[1]]$Compare.List[1:min(c(30,nrow(down.disambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
pathway.heatmap(matrixData=up.assoc,clustering = T,fontSize=5)
pathway.heatmap(matrixData=down.assoc,clustering = T,fontSize=5)