
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

•	For uniConSig analysis, the input is a target gene list provided by user containing the Entrez gene IDs. The genes included in this gene list should be functionally related to a certain disease or biological process. This gene list will be used as training gene list to calculate the functional relevance scores for all genes in the genome underlying this disease or biological process. 

•	Based on our experience, we highly recommend that the target gene list includes at least 30 genes (preferably more than 50 genes) that are highly functionally related. Random or unrelated input genes will result in ambiguous results and will make it difficult to interpret the biological meaning of the uniConSig scores. 

•	User can create a vector containing the gene list using various methods. Here is an example for creating a gene list vector containing the genes from non-homologous end joining pathway:

> myInput<-c(7518,4361,27343,27434,731751,79840,3981,2237,1791,7520,10111,2547,5591,64421)

•	To calculate uniConSig scores, please load the molecule concept dataset we compiled with precomputed Jaccard matrix:

> preCal.local<-get_data_uniConSigPreCal()

> load(preCal.local)

•	To calculate uniConSig scores for all genes in the genome, please run:

> myResult<-cal_uniConSig(myInput,preCal=preCal.data.all)

•	The object “myResult” is a “data.frame”, containing 3 columns, “gene.id”(which is entrez gene ID), “gene.name”(which is gene symbol), and “gene.score”(the uniConSig score).

•	The higher the uniConSig score of a given gene indicate the higher functional relevance of that gene underlying the function shared by the target gene list.



__E.	How to run CSEA analysis to interpret the pathways underlying a target gene list defined by genomic data__

•	To run CSEA analysis on a target experimental gene list (i.e. upregulated or downregulated genes, genetically altered genes), first build a vector for the target gene list using Entrez Gene IDs. Here is an example for building a vector for a target gene list:

> myInput<-c(3488, 1356, 10371, 9415, 8942, 718, 3157, 28983, 3638, 29113, 7035, 5507, 79971, 4250, 157570, 4233, 26830, 342667, 7042, 27346, 80339, 25984, 3776, 83439, 105374290, 2202, 23266, 51442, 2760, 51196, 57523, 121227, 65983, 121506, 727936, 57447, 2115, 10769, 84904, 4256, 254263, 8991, 283431, 501, 25840, 4602, 4998, 114907, 7368)

•	Based on our experience, we highly recommend that the target gene list includes at least 30 genes (preferably more than 50 genes).

•	Second, calculate the uniConSig scores for this target gene list:

>preCal.local<-get_data_uniConSigPreCal()

>load (preCal.local)

>myResult<-cal_uniConSig(myInput,preCal=preCal.data.all) 

•	Third, build the pathway gene sets to compare with into a list object. Here we already prebuilt two pathway gene sets into list objects, including the Hallmark and C2CP pathway gene sets from MSigdb. To run CSEA against the hallmark pathway gene sets, please run:

> data(pathway.hallmark)

> result.CSEA.hallmark<-CSEA(myResult, pathway.hallmark)

•	To run CSEA against the C2CP pathway gene sets, please run:

> data(pathway.c2cp)

> result.CSEA.c2cp<-CSEA(myResult, pathway.c2cp)

•	The result object is a data frame containing the names of the pathways included in the pathway gene sets sorted the normalized enrichment scores (NES) together with the p-values.

•	User can also define their own pathway gene sets. The object of the input pathway gene sets is a list, which can be created as below:

> pathway.my<-list()

> pathway.my[[1]]<-c(7249,3479,27330,2475,1975,1977,1978,253260,3630,5228,8408,6009,673,3091,5170,57521,207)

> pathway.my[[2]]<-c(595,596,1869,317,581,1017,1019,4193,5925,1026,7157,1647,472,7078,5111,898)

> names(pathway.my)<-c("KEGG_MTOR_SIGNALING_PATHWAY","BIOCARTA_P53_PATHWAY")

•	These codes create a list named “pathway.my”, which contains two pathways, one is consisted of entrez gene IDs of KEGG_MTOR_SIGNALING_PATHWAY from MSigdb C2CP pathway collection, and the other one BIOCARTA_P53_PATHWAY. 

•	To run CSEA analysis against user defined pathway gene sets, please run:

> result.CSEA.my<-CSEA(myResult, pathway.my)
