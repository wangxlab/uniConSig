
__Universal Concept Signature Analysis: Genome-Wide Quantification of New Biological and Pathological Functions of Genes and Pathways__

The current version is 0.99.13 (May 21st, 2019)


__A. How to install uniConSig__

1. Install the package “devtools” by running the following code in R: 
>install.packages(“devtools”)

2. You can install uniConSig from Github by running the following codes in R:
>library(devtools)
>install_github("wangxlab/uniConSig")

3. To use the uniConSig package, first load the package:
> library(uniConSig)



__B. How to run uniConSig__

For uniConSig calculation, the input is a vector containing the entrez gene IDs of the user-defined genes. User can create the vector by various methods. One simple method is to type the IDs one by one like this:

> myInput<-c(7157,4609,1:10)

The first two IDs are TP53 and MYC, the rest of IDs are just for illustration. Although there’s no restrictions on the number of genes in the input vector, we highly recommend that the input includes at least 30 genes, and the input genes should be highly related to one criterion. Random or unrelated input genes will result in the ambiguous result and it will be difficult to implicate the biological meaning of the uniConSig scores. 

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
> data(pathway.hallmark)
> result.CSEA.hallmark<-CSEA(myResult,pathway.hallmark)

The other compiled pathway collection in this package is the C2CP of MSigdb. To use C2CP pathways, run:
> data(pathway.c2cp)
> result.CSEA.c2cp<-CSEA(myResult,pathway.c2cp)

User can also define their own pathway gene sets. The object of the input pathway gene sets is a list, which can be created like this:
> pathway.my<-list()
> pathway.my[[1]]<-1:10
> pathway.my[[2]]<-5:15

These codes create a list named “pathway.my”, which contains two pathways, one is consisted of entrez gene IDs “1,2,3,4,…10”, and the other one “5,6,7,8,…15”. Then run:
> result.CSEA.my<-CSEA(myResult,pathway.my)

The result is consisted of 3 columns, the names of the pathways, the normalized enrichment scores(NES), and the p-values. You can see the first several rows of the result by;
> head(result.CSEA.my)

