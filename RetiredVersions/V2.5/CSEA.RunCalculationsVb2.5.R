##required for all calculations
source("J:/GenomeIndex/Pipeline_CSEA2/CSEA.modulesVb2.5.R")
gmtfile="J:/GenomeIndex/Pipeline_CSEA2/ConceptDb/ConceptDb20190624.rmCa.gmt"
feature.list=read_gmt(gmtfile,min=10)
feature.preCalfile=paste(gsub(".gmt","",gmtfile),".Ochiai.min10.preCal.gmt",sep="")
#batch_calSimilarity(feature.list=feature.list,feature.preCalfile=feature.preCalfile,cutoff=0.05,minsize=10,method="Ochiai")
preCalmatrix<- as.matrix(unlist(readLines(feature.preCalfile)[-1]),col=1)
###Perform CSEA for scRNAseq based on SCDE results
compare.list=c(read_gmt("J:/GenomeIndex/Pipeline_CSEA2/PathwayDb/h.all.v6.2.symbols.gmt",min=5),read_gmt("J:/GenomeIndex/Pipeline_CSEA2/PathwayDb/c2.cp.v6.2.symbols.gmt",min=5))
SCDE=read.table("J:/GenomeIndex/Pipeline_CSEA2/scDataset/scdeOut_HSCscSEQ_TrueReadCount_7447g_Quies77s_Active100s.txt",stringsAsFactors = F,fill = T,check.names = F, header = T,sep = "\t")
weight=SCDE$Signed.P.Zscore
names(weight)=SCDE$Gene
ks.result=run.weightedKS(weight,signed=F,feature.list,minsize=5,correct.overfit=FALSE,correct.outlier=FALSE,transformNegWeight=FALSE)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list,outfile,p.cut=0.01,q.cut=0.25,NES.cut=0,power=1,root=1,ECNpenalty=0.5,correct.overfit=FALSE)
up.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)
down.CSEA.result<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list,p.cut=0.05,minsize=5)
up.disambiguate<-disambiguation(CSEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,topn=min(c(30,nrow(up.CSEA.result))),p.cut=0.01)
down.disambiguate<-disambiguation(CSEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,topn=min(c(30,nrow(down.CSEA.result))),p.cut=0.01)
up.assoc <- pathwayAssociation(topPathway=up.disambiguate[[1]]$Compare.List[1:min(c(30,nrow(up.disambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
down.assoc <- pathwayAssociation(topPathway=down.disambiguate[[1]]$Compare.List[1:min(c(30,nrow(down.disambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
pathway.heatmap(matrixData=up.assoc,clustering = FALSE)
pathway.heatmap(matrixData=down.assoc,clustering = FALSE)
#################################
draw.pathway(weight=-weight,pathway="REACTOME_GENERIC_TRANSCRIPTION_PATHWAY",pathway.list=compare.list)
