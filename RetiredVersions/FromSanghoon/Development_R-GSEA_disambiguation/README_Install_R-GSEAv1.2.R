
setwd("/zfs2/xiaosongwang/sal170/18_uniConSig_CSEA2/4-1_GSEA_R_test/R-GSEA_v1.2_Test")
library("devtools")
install_github("GSEA-MSigDB/GSEA_R")

###Or 
###Optionally, a helper script has been provided to simplify use of this package. 
###Initialize the helper script with the source() command by calling 
#source(system.file('extdata', 'Run.GSEA.R', package = 'GSEA'))


## after install 
library(GSEA)
sessionInfo()

# other attached packages:
#[1] GSEA_1.2       dplyr_0.8.3    usethis_1.5.0  devtools_2.0.2

