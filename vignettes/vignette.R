## ----eval=FALSE------------------------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("uniConSigPreCal")

## ----eval=FALSE------------------------------------------------------------
#  library("uniConSigPreCal")
#  preCal.local<-get_data_uniConSigPreCal()
#  load(preCal.local)

## --------------------------------------------------------------------------
data(trList.cgc) #trList.cgc is a vector containing entrez cancer gene ids generated by Cancer Gene Census

## --------------------------------------------------------------------------
trList.my<-1:100

## ----eval=FALSE------------------------------------------------------------
#  result.cgc<-cal_uniConSig(trList.cgc,preCal=preCal.data.all)
#  result.my<-cal_uniConSig(trList.my,preCal=preCal.data.all)

## --------------------------------------------------------------------------
trList.DEG<-1:50 #user's gene list should be generated based on differentially expressed genes

## --------------------------------------------------------------------------
data(pathway.c2cp) #This package also includes a compiled list of hallmark gene sets from MSigDB

## ----eval=FALSE------------------------------------------------------------
#  result.uniConSig.DEG<-cal_uniConSig(trList.DEG,preCal=preCal.data.all)
#  result.csea.DEG.c2cp<-CSEA(result.uniConSig.DEG,pathway.c2cp) # CSEA calculation
#  head(result.csea.DEG.c2cp)

## --------------------------------------------------------------------------
sessionInfo()

