# R-GSEA 1.1 -- Gene Set Enrichment Analysis / Broad Institute 
# Currentyly, R-GSEA 1.1 was tested on R version 3.6.0

library(data.table)
library(dplyr)
library("utils")
library("tools")
library("qvalue")
library(GSEA)
#####################################################################
##### Set working directory and output directory  
setwd("R:/CaGenome/Pipelines/CSEA2/Development_R-GSEA_disambiguation")
source("R-GSEA.1.2_GeneRanking_SH2.R")
source("R-GSEA.1.2.module_disambiguation_SH2.R")
#####################################################################
##### These are for GSEA
input.ds<-"GSEAInput_ECOG_Others22s_ArmD18s_38175GeneSymb.gct"	# Input gene expression dataset file in GCT format
input.cls<-"GSEAInput_ECOG_Others22s_0_vs_ArmD18s_1.cls"		# Input class vector (phenotype) file in CLS format. Tab delimited works
gs.db<-"../PathwayDb/hallmark.c2cp.all.v6.2.GeneSymbols.gmt"		# Gene set database in GMT format
input.chip<-""   #  "GENE_SYMBOL.chip"  			#  inputchip               # CHIP File
output.directory<-"./R-GSEAout_ECOG_Others_vs_Darm"	#  Directory where to store output and results (default: "")
doc.string<-"ECOG" 
reshuffling.type      = "gene.labels"
collapse.mode         = "median"
reverse.sign          = TRUE
gsea.type  = "GSEA"
##################################################
dir.create(output.directory)


## Run GSEA: it takes about 10 minutes
GSEAout<- GSEA(   ## Input/Output Files :-------------------------------------------
	input.ds=input.ds,	# Input gene expression dataset file in GCT format
	input.cls=input.cls,		# Input class vector (phenotype) file in CLS format. Tab delimited doesn't work
	gs.db=gs.db,		# Gene set database in GMT format
	input.chip=input.chip,  			#  inputchip,               # CHIP File
	output.directory=output.directory,			# It needs "/" at the end. Directory where to store output and results (default: "")
	doc.string=doc.string,   									# Documentation string used as a prefix to name result files (default: "GSEA.analysis")
	reshuffling.type=reshuffling.type,
	collapse.mode=collapse.mode,
	reverse.sign=reverse.sign,
	gsea.type=gsea.type
)
#----------------------------------------------------------------------------------------------------------

## Run GSEA_SignalToNoise: it takes about 1~2 minutes
SignalToNoise<-GSEA_SignalToNoise(   ## Input/Output Files :-------------------------------------------
	input.ds=input.ds,	# Input gene expression dataset file in GCT format
	input.cls=input.cls	# Input class vector (phenotype) file in CLS format. Tab delimited doesn't work
)
#---------------------------

up.GSEA.result<-func_findGSEAout2(GSEAout[[1]]) # dim: 281 4   # check whether GSEAout[[2]] is factor or not.
down.GSEA.result<-func_findGSEAout2(GSEAout[[2]])

up.disambiguate<-disambiguation_RGSEA(GSEA.result=up.GSEA.result,target.score=SignalToNoise, gs.db=gs.db,upPathway=TRUE,topn=min(c(30,nrow(up.GSEA.result))),p.cut=0.01)
down.disambiguate<-disambiguation_RGSEA(GSEA.result=down.GSEA.result, target.score=SignalToNoise,gs.db=gs.db,upPathway=FALSE,topn=min(c(30,nrow(down.GSEA.result))),p.cut=0.01)
