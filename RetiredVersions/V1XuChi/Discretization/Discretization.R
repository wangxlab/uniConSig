
library(infotheo) # for discretization

######################################################################################################
setwd("C:/Users/lees18/Box Sync/lees18_UPMC_desktop/19_GeneSig_DrugResponse/5_GenoSig_Rcode/1_Discretization_ConvertToGMT_format")
inFilePath="GeneExpMatrixExample_1000g14s.txt"
outFilePath="GeneExpMatrixExample_Discretized_1000g14s.txt"
numberBins=8
######################################################################################################

data=read.table(inFilePath, sep="\t", stringsAsFactors=T, row.names=1, header=TRUE, check.names=F)

geneName <-rownames(data)
rownames(data)<-NULL
discretized=discretize(data, disc="equalwidth", nbins=numberBins)
rownames(discretized)<- geneName
write.table(discretized, outFilePath, sep="\t", row.names=T, col.names=NA, quote=F)
