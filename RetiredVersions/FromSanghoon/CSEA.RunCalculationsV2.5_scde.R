
#################################################################################
##required for all calculations
setwd("/zfs2/xiaosongwang/sal170/18_uniConSig_CSEA2/7-1_CSEAV2.3_newNorm_HSC_scSeq_UMI_nonUMI")
source("CSEA.modulesVb2.5.R")
gmtfile="../ConceptDb/ConceptDb20190624.rmCa.gmt"
readCountGCTFile<-"../../19_SymSim/4_Dataset2ObservedCounts_test/GSE68981_SCSeq_FeatureCount_HumanGeneSymb_7744g177smp.gct"
clsFile<-"../scDataset/scData_quiescent77sVsActive100s_tabdelimit.cls"
targetListFile<-"../scDataset/LimmaOutput_newNorm_scData_ObsTPMlog2_QuiescVsActive_Down500g.txt"

output.dCSEA<-"dCSEAv2.3_Output_HallmarkC2CP_newNorm_scData_TrueTPMlog2_ByLimmaDw500g.txt"

output.wCSEAUp<-"wCSEAv2.4_UpOut_HallmarkC2CP_scde_scData_TrueTPMlog2_7447g_Q77s_A100s.txt"
output.wCSEADw<-"wCSEAv2.4_DownOut_HallmarkC2CP_scde_scData_TrueTPMlog2_7447g_Q77s_A100s.txt"
###############################################################################

##required for all calculations
compare.list=c(read_gmt("../PathwayDb/h.all.v6.2.symbols.gmt",min=5),read_gmt("../PathwayDb/c2.cp.v6.2.symbols.gmt",min=5))
feature.list=read_gmt(gmtfile,min=10)
feature.preCalfile=paste(gsub(".gmt","",gmtfile),".Ochiai.min10.preCal.gmt",sep="")
preCalmatrix<- as.matrix(unlist(readLines(feature.preCalfile)[-1]),col=1)

#######################################################################
###Perform dCSEA for scRNAseq
#######################################################################

####CalUniConSig for Cancer Gene
target.data<-read.table(targetListFile,stringsAsFactors=T,header=T,sep="\t",quote="")  #15119 7
target.list<-target.data$Symbol
uniConSig=cal.uniConSig(target.list=target.list,feature.list=feature.list,preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai")
CSEA.result<-CSEA2(setNames(as.numeric(uniConSig$uniConSig), uniConSig$subjectID),compare.list,p.cut=0.05)

#write.table(CSEA.result,output.dCSEA,col.names=NA,row.names=T,sep="\t",quote=F)
#save.image("dCSEAv2.3_newNorm_scData_TrueTPMlog2_byLimmaDw500g.RData")

## PathwayAssociation
PathwayAssociationOut <- pathwayAssociation(topPathway=CSEA.result$Compare.List[1:20],compare.list,feature.list,preCalmatrix,minsize=10)
#write.table(PathwayAssociationOut ,output.up.assoc,col.names=NA,row.names=T,sep="\t",quote=F)


#######################################################################
###Perform wCSEA for scRNAseq
#######################################################################
scdeOut<-scdeDEG(readCountGCT=readCountGCTFile,clsFile=clsFile)
#write.table(scdeOut,"scdeOut_HSCscSEQ_TrueReadCount_7447g_Quies77s_Active100s.txt", col.names=NA,quote=F,sep="\t")
weight=scdeOut$Signed.P.Zscore
names(weight)=rownames(scdeOut) 

ks.result=run.weightedKS(weight,signed=F,feature.list,minsize=5,correct.overfit=FALSE,correct.outlier=FALSE,transformNegWeight=FALSE)
uniConSig.result=cal.uniConSig.ks(up.ks=ks.result[[1]],down.ks=ks.result[[2]],preCalmatrix,feature.list,outfile,p.cut=0.01,q.cut=0.25,NES.cut=0,power=1,root=1,ECNpenalty=0.5,correct.overfit=FALSE)  #correct.overfit=FALSE is default

#up.disambiguate<-disambiguation(CSEA.result=up.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=TRUE,topn=min(c(30,nrow(up.CSEA.result))),p.cut=0.01)
	down.disambiguate<-disambiguation(CSEA.result=down.CSEA.result,uniConSig.result=uniConSig.result,compare.list=compare.list,upPathways=FALSE,topn=min(c(30,nrow(down.CSEA.result))),p.cut=0.01)
	
	#up.assoc <- pathwayAssociation(topPathway=up.disambiguate[[1]]$Compare.List[1:min(c(30,nrow(up.disambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
	down.assoc <- pathwayAssociation(topPathway=down.disambiguate[[1]]$Compare.List[1:min(c(30,nrow(down.disambiguate[[1]])))],compare.list,feature.list,preCalmatrix,minsize=10)
	#pathway.heatmap(matrixData=up.assoc,clustering = FALSE)
	pathway.heatmap(matrixData=down.assoc,clustering = FALSE)
	#################################
	draw.pathway(weight=-weight,pathway="REACTOME_GENERIC_TRANSCRIPTION_PATHWAY.txt",pathway.list=compare.list)
