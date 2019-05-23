
__How to use uniConSig and CSEA R modules__

__Universal Concept Signature Analysis: Genome-Wide Quantification of New Biological and Pathological Functions of Genes and Pathways__

The current version is 0.99.13 (May 21st, 2019). The uniConSig and CSEA were built on R version 3.6.0


__A. Introduction__

•	Identifying new gene functions and pathways underlying diseases and biological processes are major challenges in genomics research. Particularly most methods for interpreting the pathways characteristic of a gene list of interest are limited by their dependence on assessing the overlapping genes or their interactome topology, which cannot account for the variety of functional relations. This is particularly problematic for pathway discovery from single cell genomics with low gene coverage, or interpreting complex pathway changes i.e. during change of cell states. 

•	Here we exploited the comprehensive sets of molecular concepts that combine ontologies, pathways, interactions, and domains to help inform the functional relations. 

•	We first developed a __universal Concept Signature (uniConSig) analysis__ for genome-wide quantification of new gene functions underlying biological or pathological processes based on the signature molecular concepts computed from known functional gene lists. 

•	We then further developed a novel __Concept Signature Enrichment Analysis (CSEA)__ for deep interpretation of the pathways enriched in a gene list of interest defined by genomic data. This method is grounded on the framework of shared concept signatures between gene sets at multiple functional levels, thus overcomes the limitations of the current methods. 

•	Through meta-analysis of transcriptomic datasets of cancer cell line models and single hematopoietic stem cells, we demonstrate the broad applications of CSEA on pathway discovery from gene expression and single cell transcriptomic datasets for genetic perturbations and change of cell states.


__B.	Basic requirements__

•	The uniConSig and CSEA modules are compiled in an R package “uniConSig” held at https://github.com/wangxlab/uniConSig. 
•	To install the package, first install R from CRAN: https://cran.r-project.org/
•	For more user-friendly interface, R-Studio can be installed from here: https://www.rstudio.com/products/rstudio/download/


__C. How to install uniConSig__

•	First, please install the package “devtools” by running the following code in R:
>install.packages(“devtools”)

•	Second, please install uniConSig from Github by running the following codes in R:
>library(devtools)<br />
>install_github("wangxlab/uniConSig")

•	To use the uniConSig package, first load the package:
> library(uniConSig)


__D.	How to run uniConSig analysis for the discovery of new gene functions underlying certain diseases or biological processes__

For uniConSig calculation, the input is a vector containing the entrez gene IDs of the user-defined genes. User can create the vector by various methods. One simple method is to type the IDs one by one like this:

> myInput<-c(7518,4361,27343,27434,731751,79840,3981,2237,1791,7520,10111,2547,5591,64421)

This gene set was taken from MSigdb’s pathway collection “KEGG_NON_HOMOLOGOUS_END_JOINING”. Although there’s no restrictions on the number of genes in the input vector, we highly recommend that the input includes at least 30 genes, and the input genes should be highly related to one criterion. Random or unrelated input genes will result in the ambiguous result and it will be difficult to implicate the biological meaning of the uniConSig scores. 

To calculate uniConSig scores, simply run:
> myResult<-cal_uniConSig(myInput)

The “myResult” will be a “data.frame”, containing 3 columns, “gene.id”(which is entrez gene ID), “gene.name”(which is gene symbol), and “gene.score”(the uniConSig score). Note that this is only for illustration purpose, so the data set used for this calculation is only a small part of the total data set. 

To obtain the total dataset, run:
> preCal.local<-get_data_uniConSigPreCal()

Then load the data into environment:
> load(preCal.local)

To run uniConSig based on the total dataset, run:
> myResult<-cal_uniConSig(myInput,preCal=preCal.data.all)




__C. How to run CSEA__

The result of uniConSig can be directly put into CSEA calculation. Before calculation of CSEA, we have to define the pathways(gene sets) that we want to test. To use the hallmark gene sets from MSigdb, for example, run:
> data(pathway.hallmark)<br />
> result.CSEA.hallmark<-CSEA(myResult,pathway.hallmark)<br />

The other compiled pathway collection in this package is the C2CP of MSigdb. To use C2CP pathways, run:<br />
> data(pathway.c2cp)<br />
> result.CSEA.c2cp<-CSEA(myResult,pathway.c2cp)<br />

User can also define their own pathway gene sets. The object of the input pathway gene sets is a list, which can be created like this:
> pathway.my<-list()<br />
> pathway.my[[1]]<-c(7249,3479,27330,2475,1975,1977,1978,253260,3630,5228,8408,6009,673,3091,5170,57521,207)<br />
> pathway.my[[2]]<-c(595,596,1869,317,581,1017,1019,4193,5925,1026,7157,1647,472,7078,5111,898)<br />
> names(pathway.my)<-c("KEGG_MTOR_SIGNALING_PATHWAY","BIOCARTA_P53_PATHWAY")<br />

These codes create a list named “pathway.my”, which contains two pathways, one is consisted of entrez gene IDs KEGG_MTOR_SIGNALING_PATHWAY from MSigdb C2CP pathway collection, and the other one BIOCARTA_P53_PATHWAY“5,6,7,8,…15”. Then run:
> result.CSEA.my<-CSEA(myResult,pathway.my)

[1] "Calculating permutations..."<br />
[1] "Done."<br />
[1] "Calculating NES"<br />


The result is consisted of 3 columns, the names of the pathways, the normalized enrichment scores(NES), and the p-values. You can see the first several rows of the result by;
> head(result.CSEA.my)
                   
                        names      NES  pValue   
1 BIOCARTA_P53_PATHWAY 1.932083  0.000 <br />
2 KEGG_MTOR_SIGNALING_PATHWAY 1.357480  0.057 <br />


