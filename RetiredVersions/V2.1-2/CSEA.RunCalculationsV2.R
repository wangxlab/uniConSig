##required for all calculations
source("Z:/CaGenome/Pipelines/CSEA2/CSEA.modulesV2.1.R")
gmtfile="Z:/CaGenome/Pipelines/CSEA2/ConceptDb/ConceptDb20190624.rmCa.gmt"
feature.list=read_gmt(gmtfile,min=10)
feature.preCalfile=paste(gsub(".gmt","",gmtfile),".Ochiai.min10.preCal.gmt",sep="")
#batch_calSimilarity(feature.list=feature.list,feature.preCalfile=feature.preCalfile,cutoff=0.05,minsize=10,method="Ochiai")
preCalmatrix<- as.matrix(unlist(readLines(feature.preCalfile)[-1]),col=1)

###Perform D-CSEA2 for dichotomous gene list from scDataset
compare.list=c(read_gmt("Z:/CaGenome/Pipelines/CSEA2/PathwayDb/h.all.v6.2.symbols.gmt",min=5),read_gmt("J:/GenomeIndex/Pipeline_CSEA2/PathwayDb/c2.cp.v6.2.symbols.gmt",min=5))
target.list<-read.table ("Z:/CaGenome/Pipelines/CSEA2/scDataset/downGenes_HSC_hgncSym.txt",header=FALSE,sep="\t")$V1
uniConSig=cal.uniConSig(target.list=target.list,feature.list=feature.list,preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai")
CSEA.result<-CSEA2(setNames(as.numeric(uniConSig$uniConSig), uniConSig$subjectID),compare.list,p.cut=0.05)

###Perform W-CSEA2 for scRNAseq
compare.list=c(read_gmt("Z:/CaGenome/Pipelines/CSEA2/PathwayDb/h.all.v6.2.symbols.gmt",min=5),read_gmt("J:/GenomeIndex/Pipeline_CSEA2/PathwayDb/c2.cp.v6.2.symbols.gmt",min=5))
limma=limmaDGE(gctFile="Z:/CaGenome/Pipelines/CSEA2/scDataset/scData_quiescentVsActive.gct",clsFile="Z:/CaGenome/Pipelines/CSEA2/scDataset/scData_quiescentVsActive.cls")
weight=limma$T.Value
names(weight)=row.names(limma)
ks.result=run.weightedKS(weight,feature.list,minsize=10,correct.overfit=FALSE,correct.outlier=FALSE)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list,outfile,p.cut=0.01,q.cut=0.25,NES.cut=0,power=1,root=1,ECNpenalty=0.5,correct.overfit=FALSE)
up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)
