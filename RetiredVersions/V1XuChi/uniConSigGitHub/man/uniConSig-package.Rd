\name{uniConSig-package}
\alias{uniConSig-package}
\alias{uniConSig}
\docType{package}
\title{
Calculate uniConSig scores using pre-calculated files
}
\description{
\packageDescription{uniConSig}
}
\details{
The package include two major functions:

\code{cal_uniConSig(trList,preCal=preCal.data)}

\code{CSEA(result.uniConSig,posList,nPermu=1000)}

The input of uniConSig is trList, which is a vector containing the entrez gene IDs of interest. The output of cal_uniConSig can be directly put into the second function "CSEA" as the first argument. The second and the third argument specify the pathway gene set. Currently only c2cp and hallmark pathways from MSigdb are compiled.
}
\author{
\packageAuthor{uniConSig}

Maintainer: \packageMaintainer{uniConSig}
}
\references{
https://github.com/wangxlab/uniConSig
}
\keyword{ package }
\examples{
#See vignettes for full calculation of uniConSig scores.
trList.call<-1:100
result<-cal_uniConSig(trList.call)
pathway.my<-list()
pathway.my[[1]]<-1:10
pathway.my[[2]]<-5:15
result.csea<-CSEA(result,pathway.my)
}
