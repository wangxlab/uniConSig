##required for all calculations
source("R:/CaGenome/Pipelines/CSEA2/CSEA.modulesV2.2.R")
gmtfile="R:/CaGenome/Pipelines/CSEA2/ConceptDb/ConceptDb20190624.rmCa.gmt"
feature.list=read_gmt(gmtfile,min=10)
feature.preCalfile=paste(gsub(".gmt","",gmtfile),".Ochiai.min10.preCal.gmt",sep="")
#batch_calSimilarity(feature.list=feature.list,feature.preCalfile=feature.preCalfile,cutoff=0.05,minsize=10,method="Ochiai")
preCalmatrix<- as.matrix(unlist(readLines(feature.preCalfile)[-1]),col=1)


####CalUniConSig for Cancer Gene
CGC=read.table("R:/CaGenome/Pipelines/CSEA2/CancerGene/cancer_gene_census.csv",header=TRUE,sep=",")
ca.uniConSig.all=cal.uniConSig(target.list=CGC$Gene.Symbol,feature.list,preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai") #method=="Ochiai" or "Jaccard"
benchmark=foldCrossVal.benchmarkUniConSig(CGC$Gene.Symbol,feaure.list,preCalmatrix,fold=5,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai")

###Perform CSEA for dichotomous gene list from scDataset
compare.list=c(read_gmt("R:/CaGenome/Pipelines/CSEA2/PathwayDb/h.all.v6.2.symbols.gmt",min=10),read_gmt("R:/CaGenome/Pipelines/CSEA2/PathwayDb/c2.cp.v6.2.symbols.gmt",min=10))
target.list<-read.table ("R:/CaGenome/Pipelines/CSEA2/scDataset/downGenes_HSC_hgncSym.txt",header=FALSE,sep="\t")$V1
uniConSig=cal.uniConSig(target.list=target.list,feature.list=feature.list,preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai")
CSEA.result<-CSEA2(setNames(as.numeric(uniConSig$uniConSig), uniConSig$subjectID),compare.list,p.cut=0.05)

###Perform CSEA for sorted gene list
compare.list=c(read_gmt("R:/CaGenome/Pipelines/CSEA2/PathwayDb/h.all.v6.2.symbols.gmt",min=10),read_gmt("R:/CaGenome/Pipelines/CSEA2/PathwayDb/c2.cp.v6.2.symbols.gmt",min=10))
limma=read.table("R:/CaGenome/Pipelines/CSEA2/CaKDDataset/DifferentialExpressionAnalysisByLimma_MDA-468_TP53KO.tsv",header=TRUE,sep="\t")
weight=limma$T.Value
names(weight)=limma$Gene.Sym
ks.result=run.weightedKS(weight,feature.list,minsize=10,correct.overfit=FALSE)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list,outfile,p.cut=0.01,q.cut=0.25,NES.cut=0,power=1,root=1,ECNpenalty=0.5,correct.overfit=FALSE)
up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05)
up.disambiguate<-disambiguation(CSEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,topn=30,p.cut=0.01)
down.disambiguate<-disambiguation(CSEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,topn=30,p.cut=0.01)

###Perform CSEA for scRNAseq
compare.list=c(read_gmt("R:/CaGenome/Pipelines/CSEA2/PathwayDb/h.all.v6.2.symbols.gmt",min=5),read_gmt("R:/CaGenome/Pipelines/CSEA2/PathwayDb/c2.cp.v6.2.symbols.gmt",min=5))
limma=limmaDGE(gctFile="R:/CaGenome/Pipelines/CSEA2/scDataset/scData_quiescentVsActive.gct",clsFile="R:/CaGenome/Pipelines/CSEA2/scDataset/scData_quiescentVsActive.cls")
weight=limma$T.Value
names(weight)=row.names(limma)
ks.result=run.weightedKS(weight,feature.list,minsize=10,correct.overfit=FALSE,correct.outlier=TRUE)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list,outfile,p.cut=0.01,q.cut=0.25,NES.cut=0,power=1,root=1,ECNpenalty=0.5,correct.overfit=FALSE)
up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)
up.disambiguate<-disambiguation(CSEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,topn=30,p.cut=0.01)
down.disambiguate<-disambiguation(CSEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,topn=30,p.cut=0.01)
up.assoc <- pathwayAssociation(topPathway=up.CSEA.result$Compare.List[1:30],compare.list,feature.list,preCalmatrix,minsize=10)
down.assoc <- pathwayAssociation(topPathway=down.CSEA.result$Compare.List[1:20],compare.list,feature.list,preCalmatrix,minsize=10) # the top pathway number should be less than the total number of significant pathways
