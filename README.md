# How to use uniConSig and CSEA R modules

__Universal Concept Signature Analysis: Genome-Wide Quantification of New Biological and Pathological Functions of Genes and Pathways__

The current version is 2.0

•	NOTE: The methods for CSEA2.0 are developed in our most recent study which have been incorporated into a new "indepthPathway" package: https://github.com/wangxlab/indepthPathway

•	Please CITE both of the following publications:

•	1. Xu Chi, Maureen A. Sartor, Sanghoon Lee, Meenakshi Anurag, Snehal Patil, Pelle Hall, Matthew Wexler, Xiaosong Wang#. Universal Concept Signature Analysis: Genome-Wide Quantification of New Biological and Pathological Functions of Genes and Pathways. Briefings in Bioinformatics, 2020 Sep 25;21(5):1717-1732.

•	2. Lee S, Deng L, Wang Y, Wang K, Sartor MA, Wang XS# IndepthPathway: an integrated tool for in-depth pathway enrichment analysis based on single-cell sequencing data. Bioinformatics. 2023 Jun 1;39(6): btad325.

## A. Introduction

•	Identifying new gene functions and pathways underlying diseases and biological processes are major challenges in genomics research. Particularly most methods for interpreting the pathways characteristic of a gene list of interest are limited by their dependence on assessing the overlapping genes or their interactome topology, which cannot account for the variety of functional relations. This is particularly problematic for pathway discovery from single cell genomics with low gene coverage, or interpreting complex pathway changes i.e. during change of cell states. 

•	Here we exploited the comprehensive sets of molecular concepts that combine ontologies, pathways, interactions, and domains to help inform the functional relations. 

•	We first developed a __universal Concept Signature (uniConSig) analysis__ for genome-wide quantification of new gene functions underlying biological or pathological processes based on the signature molecular concepts computed from known functional gene lists. 

•	We then further developed a novel __Concept Signature Enrichment Analysis (CSEA)__ for deep interpretation of the pathways enriched in a gene list of interest defined by genomic data. This method is grounded on the framework of shared concept signatures between gene sets at multiple functional levels, thus overcomes the limitations of the current methods. 

•	Through meta-analysis of transcriptomic datasets of cancer cell line models and single hematopoietic stem cells, we demonstrate the broad applications of CSEA on pathway discovery from gene expression and single cell transcriptomic datasets for genetic perturbations and change of cell states.

## B.	Basic requirements

•	The uniConSig and CSEA modules are compiled in an R package “uniConSig” held at https://github.com/wangxlab/uniConSig. 

•	To install the package, first install R from CRAN: https://cran.r-project.org/

•	For more user-friendly interface, R-Studio can be installed from here: https://www.rstudio.com/products/rstudio/download/


## C. Installing R Packages

•	Find "uniConSig_CSEA.modulesV2.0.R" file in your folder and see line 1. There is the list of R packages you may have to install unless you did.
```
packages=c("stringr","qvalue","vegan","tidyr","moments","limma","dplyr","gplots","RColorBrewer","corrplot","pheatmap","igraph","otuSummary") 
sapply(packages,require,character=TRUE)
```

## D.	How to calculate uniConSig score and run D-CSEA for pathway enrichment analysis
•	To calculate uniConSig scores, please load the molecule concept dataset we compiled with precomputed Jaccard matrix:
```
#IMPORTANT: this code block is required for all calculations including both uniConSig and CSEA.
setwd("your working directory")
source("uniConSig_CSEA.modulesVb2.8.R")
gmtfile="ConceptDb/ConceptDb20190624.rmCa.gmt"
feature.list=read_concepts(gmtfile)

#if you are using your own molecular concept database, please use the following code generate preCal file containing information about concept redundancy. The molecule concepts need to be provided as gmt file format.
#batch_calSimilarity(feature.list=feature.list,feature.preCalfile=feature.preCalfile)
#Load precomputed concept redundancy data generated from the above code
feature.preCalfile=paste(gsub(".gmt","",gmtfile),".Ochiai.min10.preCal.gmt",sep="")
preCalmatrix<- as.matrix(unlist(readLines(feature.preCalfile)[-1]),col=1)
```

## Part I. Calculate uniConSig scores to compute a functional relevance of human genes underlying a target gene list such as known cancer genes 
```
#rm.overfit should be set to TRUE for calculating gene set scores

####CalUniConSig for Cancer Gene
CGC=read.table("GeneCensus/CancerGeneCensus_allSat Jan 30 17 31 29 2021.tsv",header=TRUE,sep="\t")
ca.uniConSig.all=cal.uniConSig.cut(target.list=CGC$Gene.Symbol,feature.list,preCalmatrix,rm.overfit=TRUE)
write.table(ca.uniConSig.all,file="GeneCensus/Cancer.UniConSig2.1.PANCAN.tsv",quote=F,row.names=FALSE,col.names = TRUE,sep="\t")
#optional: benchmark cancer gene prioritization results
benchmark=foldCrossVal.benchmarkUniConSig(CGC$Gene.Symbol,feaure.list,preCalmatrix,fold=5)

####CalUniConSig for Immune Gene
immune=read.table("GeneCensus/IMMUNE_SYSTEM_PROCESS.tsv",skip=2,sep="\t")
im.uniConSig.all=cal.uniConSig.cut(target.list=immune$V1,feature.list,preCalmatrix,rm.overfit=TRUE)
write.table(im.uniConSig.all,file="GeneCensus/Immune.UniConSig2.1.tsv",quote=F,row.names=FALSE,col.names = TRUE,sep="\t")

```
## Part II. perform D-CSEA for pathway enrichment analysis

```
###Perform D-CSEA for dichotomous experimental gene list from scDataset
compare.list=c(read_concepts("PathwayDb/h.all.v7.5.1.symbols.gmt"),read_concepts("PathwayDb/c2.cp.v7.5.1.symbols.gmt"))
target.list<-read.table ("scDataset/downGenes_HSC_hgncSym.txt",header=FALSE,sep="\t")$V1
#perform deep functional interpretation of the target gene list and calculate uniConSig scores.The parameter rm.overfit should set as false for pathway enrichment analysis, which will give high weights for the genes included in the experimental gene list
uniConSig=cal.uniConSig(target.list=target.list,feature.list=feature.list,preCalmatrix,rm.overfit=F)
CSEA.result<-CSEA2(setNames(as.numeric(uniConSig$uniConSig), uniConSig$subjectID),compare.list,p.cut=0.05)#p.cut: the p value cutoff for significant pathways
#disambiguate top enriched pathways.
topn=100 #specify the number of top pathways to disambiguate
disambiguate<-disambiguation.CSEA(GEA.result=CSEA.result,uniConSig.result=uniConSig,compare.list=compare.list,topn=min(c(topn,nrow(CSEA.result))),p.cut=0.01) ##upPathways use TRUE/FALSE for W-CSEA or use NULL D-CSEA
#compute functional associations between selected top pathways
selectn=30 #specify the number of top pathways to compute associations
assoc <- pathwayAssociation(topPathway=disambiguate[[1]]$Compare.List[1:min(c(selectn,nrow(disambiguate[[1]])))],compare.list,feature.list,preCalmatrix)
#draw heatmaps showing the functional associations between selected top pathways 
pdf(file="D-CSEA-PathwayAssocHeatmapNetwork_20220627.pdf",width=10, height=10)
pathway.heatmap(matrixData=assoc,clustering = TRUE,fontSize=8)
#draw network figure. NES.cut determines the levels of significance for the pathway similarity to be shown as edges. The higher the NES, the less connections will be shown in the network.
draw.network(pathway.out=disambiguate[[1]],assoc=assoc,NES.cut=2)
graphics.off()
```
