library(limma)
limmaDGE<-function(gctFile,clsFile){ # the comparison is based on increasingly sorted orders of group labels with positive t value showing upreulated in class 1
  ############################################
  ## Gene expression matrix data process 
  ############################################
  myTableDF=read.table(gctFile, sep="\t", stringsAsFactors=F, header=F, row.names=NULL, check.names=F, fill=TRUE)[-1,-2]
  myTableNoDup<-myTableDF[which(!duplicated(myTableDF[,1])),]#REVISE
  colnames(myTableNoDup)<-myTableNoDup[1,]
  rownames(myTableNoDup)<-myTableNoDup[,1]
  ## convert character to numeric in data.frame
  myTableNumericDF<-as.data.frame(sapply(myTableNoDup[-1,-1], as.numeric))
  rownames(myTableNumericDF)<-rownames(myTableNoDup[-1,-1])
  ## Remove genes (rows) that have zero variance
  myTableNumDF_noZeroGene<-myTableNumericDF[(which(!apply(myTableNumericDF,1,var)==0)),]
  
  ############################################
  ## Get the class factors.
  ############################################
  myClassFactorDF=read.table(clsFile, sep="\t", stringsAsFactors=F, header=F, row.names=NULL, check.names=F, fill=TRUE)[-1,]
  fac<-as.factor(myClassFactorDF[1,])
  label=sort(unique(as.character(fac)))
  ###########################################
  ## Calculate log2 fold change
  ###########################################
  log2Foldchange<-apply(myTableNumDF_noZeroGene,1,function(x){
  log2(mean(as.numeric(x[which(fac==label[2])]))/mean(as.numeric(x[which(fac==label[1])])))
  })
  log2FC_DF<-as.data.frame(matrix(log2Foldchange,ncol=1))
  rownames(log2FC_DF)<-names(log2Foldchange)
  colnames(log2FC_DF)<-"log2FC"
  
  ###########################################
  ## Run limma
  ###########################################
  fit<-lmFit(myTableNumDF_noZeroGene,design=model.matrix(~fac))
  fit.eBayes <- eBayes(fit)
  limmaOut<-topTable(fit.eBayes,number=100000000, adjust.method="BH")  # BH, Benjamini and Hochberg FDR adjusted p-value.
  
  ###########################################
  ## Merge Log2FC and Limma output data
  ###########################################
  limmaOutLog2FC<-merge(log2FC_DF,limmaOut, by=0)
  rownames(limmaOutLog2FC)<-limmaOutLog2FC[,1]
  limmaOutLog2FC$Row.names<-NULL
  limmaOutLog2FC$logFC<-NULL
  limmaOutLog2FC$B<-NULL
  colnames(limmaOutLog2FC)<-c("Log2FC","Average Expression", "T.Value","P.Value","FDR_Q.Value")
  limmaOutLog2FC[limmaOutLog2FC=="-Inf"|limmaOutLog2FC=="Inf"]<-NA
  limmaOutLog2FCsorted<-limmaOutLog2FC[order(limmaOutLog2FC$FDR_Q.Value),]
  return(limmaOutLog2FCsorted)
}