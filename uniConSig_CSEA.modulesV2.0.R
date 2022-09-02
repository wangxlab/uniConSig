packages=c("stringr","qvalue","vegan","tidyr","moments","limma","dplyr","gplots",
           "RColorBrewer","corrplot","pheatmap","igraph","otuSummary","pROC",
           "matrixStats","pacman") 
sapply(packages,require,character=TRUE)

library(stringr)
library(vegan)
library(qvalue)
library(tidyr)
library(moments)
library(limma)
library(dplyr)
library(gplots)
library(RColorBrewer)
library(corrplot)
library(igraph)
library(otuSummary)
normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
#filter by one list
filter <- function (x,filter,inlist) {
  if (inlist){
    x <- x[x %in% filter]  
  }
  else{
    x <- x[!x %in% filter]  
  }
  return (x)
}

#filter by two list
filter.double <- function (x,include,exclude) {
  x <- x[(x %in% include) & (!x %in% exclude)]  
  return (x)
}

#turn genotype list into a matrix fast
genolist2matrix<-function(genolist){
  un.genolist <- unlist(genolist)
  cell<-sort(unique(un.genolist))
  tmp.bin <- matrix(0, nrow = length(genolist), ncol = length(cell))
  dimnames(tmp.bin) <- list(names(genolist), cell)
  ij <- cbind(rep(1:length(genolist), lengths(genolist)), match(un.genolist,cell))
  tmp.bin[ij] <- 1
  return(tmp.bin)
}
#turn matrix into genolist file
matrix2genofile<-function(matrix,outfile){
  for (i in 1:nrow(matrix)){
    tmp.out=c(unlist(strsplit(rownames(matrix)[i],":")),names(which(matrix[i,]==1)))
    write(tmp.out,file=outfile,append=TRUE,sep="\t",ncolumns = length(tmp.out))
  }
}

#take nth root
nthroot<-function(x,n){
  abs(x)^(1/n)*sign(x)
}


read_gmt <- function(fileName, min = 5) {
  Lines <- fread(fileName, sep = "")[[1]] #read lines. This seems to remove the annotation line start with #
  read.list <- lapply (Lines, function(x) {
    genotype=strsplit(x, "\t")[[1]]
    gname=paste(genotype[1],genotype[2],sep=":")
    return(list(gname,genotype[-c(1:2)]))
  })
  genotype.list=lapply(read.list, `[[`, 2) 
  names(genotype.list)= lapply(read.list, `[[`, 1)
  genotype.list=genotype.list[lengths(genotype.list)>=min]
  return(genotype.list)
}
write_gmt <- function(genotype.list,type,file.name) {
  if (file.exists(file.name)) {
    message(paste(file.name, "exists, deleting..."))
    file.remove(file.name)
  }
  write("#GMT file written from a genotypelist", file=file.name, sep="\t",append=TRUE, ncolumns=1)
  for (i in seq_along(genotype.list)) {
    genotype <- genotype.list[[i]]
    gname=names(genotype.list)[i]
    ncolumns <- 2 + length(genotype)
    write(c(gname,type,genotype), file=file.name, sep="\t",append=TRUE, ncolumns=ncolumns)
  }
}

require(data.table)
read_concepts <- function(fileName, min = 10) { # min: set the cutoff for minimal number of genes in the concepts or pathways
  Lines <- fread(fileName, sep = "")[[1]] #read lines. This seems to remove the annotation line start with #
  read.list <- lapply (Lines, function(x) {
    genotype=strsplit(x, "\t")[[1]]
    gname=genotype[1]
    return(list(gname,genotype[-c(1:2)]))
  })
  genotype.list=lapply(read.list, `[[`, 2) 
  names(genotype.list)= lapply(read.list, `[[`, 1)
  genotype.list=genotype.list[lengths(genotype.list)>=min]
  return(genotype.list)
}
#find concepts that contains the cosmicID in a conceptList
findConcepts<-function(cosmicID,conceptList){
  concepts<-c()
  for(i in 1:length(conceptList)){
    if(cosmicID %in% conceptList[[i]]){
      concepts<-append(concepts,names(conceptList)[i])
    }
  }
  return(concepts)
}
# generate fold cross validation groups for a target list.

foldCrossGroups<-function(target.list,fold=5){
  folds <- cut(seq(1,length(target.list)),breaks=fold,labels=FALSE)
  testSet=list()
  for(i in 1:5){
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testSet[i]<-list(i=target.list[testIndexes])
  } 
  return(testSet)
}

#fast calculate feature redundancy
calSimilarity<-function(subject.call,feature.list,cutoff=0.1,method="Ochiai"){ #method=="Ochiai" or "Jaccard"
  tmp.feature<-findConcepts(subject.call,feature.list)
  tmp.feature<-tmp.feature[tmp.feature %in% names(feature.list)]
  tmp.featurelist<-feature.list[tmp.feature]
  if (length(tmp.featurelist)==0){
    tmp.out<-c(as.character(subject.call),0)
  }else{
    tmp.featurelist.bin<-genolist2matrix(tmp.featurelist)
    if (method=="Ochiai"){
      sim.matrix<-1-as.matrix(designdist(tmp.featurelist.bin, method = "1-J/sqrt(A*B)"))
    }else if(method=="Jaccard"){
      sim.matrix<-1-as.matrix(designdist(tmp.featurelist.bin, method = "(A+B-2*J)/(A+B-J)"))
    }else if(method=="Raup-Crick"){
      sim.matrix<-1-as.matrix(designdist(tmp.featurelist.bin, method = "1-phyper(J-1, A, P-A, B)"))
    }else if (method=="Kulczynski"){
      sim.matrix<-1-as.matrix(designdist(tmp.featurelist.bin, method = "1-(J/2)*(1/A+1/B)"))
    }
    sim.matrix <- ifelse(sim.matrix<cutoff,0,sim.matrix)
    sum.similarity<-colSums(sim.matrix)
    ECN<-sum(1/sqrt(sum.similarity))
    tmp.out<-c(as.character(subject.call),ECN)
    tmp.out<-append(tmp.out,paste(names(sum.similarity),sum.similarity,sep="@"))
  }
  return(tmp.out)
}

#precalculate feature redundancy. Default parameters: cutoff=0.05,minsize=10,method="Ochiai"
batch_calSimilarity<-function(feature.list,feature.preCalfile,cutoff=0.05,minsize=10,method="Ochiai"){
  feature.list<-feature.list[which(sapply(feature.list,function(x) length(x)>=minsize))]
  subject.id<-unique(unlist(feature.list))
  suppressWarnings(file.remove(feature.preCalfile))
  cat(paste("#Precalculations of feature redundancies based on ",method," Similarity Matrix with ",method," Index cut off of ",cutoff,", min feature size of ",minsize,"\n",sep=""),file=feature.preCalfile,append=TRUE)
  for (i in 1:length(subject.id)){
    tmp.out<-calSimilarity(subject.call=subject.id[i],feature.list=feature.list,cutoff=cutoff,method=method)
    write(tmp.out,file=feature.preCalfile,append=TRUE,sep="\t",ncolumns = length(tmp.out))
    if(i %% 50==0){
      print(paste("Processed ",i," subject ids",sep=""))
    }
  }
}

#Calculate weight based on Ochiai Index
CalWeight<-function(list1,list2,method="Ochiai"){ #method=="Ochiai" or "Jaccard"
  tmp.intersect<-intersect(list1,list2)
  tmp.union<-union(list1,list2)
  if (method=="Ochiai"){
    tmp.weight=length(tmp.intersect)/sqrt(length(list1)*length(list2))
    tmp.weight.1=(length(tmp.intersect)-1)/sqrt(length(list1)*length(list2)) 
  }else if (method=="Jaccard"){
    tmp.weight=length(tmp.intersect)/length(tmp.union)
    tmp.weight.1=(length(tmp.intersect)-1)/length(tmp.union)
  }
  return(c(tmp.weight,tmp.weight.1,length(tmp.intersect),length(list1),length(list2)))
}


#batch calculate weight for a target list and a compendia of comparing lists
batch_CalWeight<-function(target.list,compare.list,method="Ochiai"){ #method=="Ochiai" or "Jaccard"
  target.list=as.character(target.list)
  target.result<-data.frame(matrix(ncol = 6, nrow = 0))
  for (i in 1:length(compare.list)){
    target.result[nrow(target.result) + 1,]<-c(names(compare.list)[i],CalWeight(target.list,compare.list[[i]],method=method))
  }
  colnames(target.result)<-c("Feature.List","Weight","Weight.1","Intersect","Target.Size","Compare.Size")
  return(target.result)
}

#calculate uniConSig scores. Default parameters: Selecting Concepts: minsize=10,weight.cut=0.05; uniConSig algorithm: power=1,root=1,ECNpenalty=0.5,method="Ochiai"
cal.uniConSig<-function(target.list,feature.list,preCalmatrix,minsize=10,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai",rm.overfit){ #method=="Ochiai" or "Jaccard"
  if (!exists("feature.list")){
    stop("please provide the list of concepts to variable: feature.list")
  }else if (!exists("preCalmatrix")){
    stop("please provide the preCalmatrix for concepts to variable: preCalmatrix")
  }
  target.weight=batch_CalWeight(target.list=target.list,compare.list=feature.list,method=method)
  target.weight=target.weight[as.numeric(target.weight$Compare.Size)>minsize&as.numeric(target.weight$Weight)>weight.cut,]
  if (nrow(target.weight)<10){
    print(paste("There are only",nrow(target.weight),"signature features with weights, this means that the target list is functionally heterogeneous, please lower the cutoff for weight.cut"))
    return(NULL)
  }else{
    result<-data.frame(matrix(ncol = 3, nrow = 0))
    for (i in 1:nrow(preCalmatrix)){
      tmp.line<-unlist(strsplit(preCalmatrix[i,],"\t"))
      subjectID=as.character(tmp.line[1])
      if(length(tmp.line)==2){
        result[nrow(result) + 1,]<-c(tmp.line[1],0)
      }
      tmp.strsplit=strsplit(tmp.line[3:length(tmp.line)],"@")
      tmp.strsplit=tmp.strsplit[which(lengths(tmp.strsplit)==2)]
      tmp.epsilon<-as.data.frame(do.call(rbind, tmp.strsplit),stringsAsFactors=FALSE)
      colnames(tmp.epsilon)<-c("Feature.List","Epsilon")
      tmp.epsilon[,"Epsilon"]=(as.numeric(tmp.epsilon[,"Epsilon"]))^root
      
      #number of concepts/ mean of epsilon
      ECN=nrow(tmp.epsilon)/mean(as.numeric(tmp.epsilon[,"Epsilon"]),trim=0.3)
      tmp.data<-merge(tmp.epsilon,target.weight,by.x="Feature.List",by.y="Feature.List")
      tmp.data$Epsilon=as.character(tmp.data$Epsilon)
      tmp.data$Feature.List=as.character(tmp.data$Feature.List)
      if (nrow(tmp.data)==0){
        next
      }
      if (rm.overfit==TRUE){
        if (subjectID %in% target.list){
          tmp.data$Weight4Cal=tmp.data$Weight.1
        }else{
          tmp.data$Weight4Cal=tmp.data$Weight
        }
      }else{
        tmp.data$Weight4Cal=tmp.data$Weight
      }
      tmp.data$Epsilon=as.numeric(as.character(tmp.data$Epsilon))
      tmp.data$Weight4Cal=as.numeric(as.character(tmp.data$Weight4Cal))
      #ECN=sum(1/(tmp.data$Epsilon^root))
      uniConSig<-sum((tmp.data$Weight4Cal^power)/(tmp.data$Epsilon))/(ECN^ECNpenalty)
      result[nrow(result) + 1,]<-c(tmp.line[1],uniConSig,ifelse(subjectID %in% target.list,1,0))
      if(i %% 5000==0){
        print(paste("Processed ",i," subject IDs",sep=""))
      }
    }
    colnames(result)<-c("subjectID","uniConSig","Target.List")
    result$uniConSig=normalize(as.numeric(as.character(result$uniConSig)))
    result=result[order(result$uniConSig, decreasing = TRUE),]
    rownames(result)=1:nrow(result)
    return(result)   
  }
}
library(pROC)
cal.uniConSig.cut<-function(target.list,feature.list,preCalmatrix,minsize=10,rm.overfit=TRUE){
  uniConSig.all=cal.uniConSig(target.list=target.list,feature.list=feature.list,preCalmatrix=preCalmatrix,minsize=minsize,rm.overfit=rm.overfit) #method=="Ochiai" or "Jaccard"
  test.roc<-roc(uniConSig.all$Target.List,uniConSig.all$uniConSig,levels=c(0,1),quiet=FALSE,direction = "<")
  opt.cut<-coords(test.roc, x="best", input="threshold", best.method="youden",best.weights = c(1, 0.5),transpose=TRUE)
  uniConSig.all$opt.cut=rep(opt.cut["threshold"],nrow(uniConSig.all))
  colnames(uniConSig.all)[1]="GeneSym"
  return(uniConSig.all)
}
#calculated weighted KS using weightedKSV2 based on a gene list with t-values as weights (such as the gene list sorted by t-values)
weightedKS.batch<-function(weight,feature.list,preCal=TRUE,high.enrich=TRUE,p.cut=0.05,correct.overfit=FALSE,correct.outlier=FALSE,minsize=5){
  feature.list<-lapply(feature.list,filter,filter=names(weight),inlist=TRUE)
  feature.list<-feature.list[which(sapply(feature.list,function(x) length(x)>=minsize))]
  if(high.enrich){
    print("Using Weight to rank subjects")
    tmp.weight<-sort(weight) ##use 1-AUC to favor the sensitive cell lines
  }else{
    print("Using 1-weight to rank subjects")
    tmp.weight=sort(max(weight)-weight)
  }
  if (correct.outlier){
    tmp.weight=tmp.weight-min(tmp.weight)
    mySlope=max(tmp.weight)/which.max(tmp.weight)
    max.index=which.max((mySlope*1:length(tmp.weight)-tmp.weight)/sqrt(1+mySlope^2))
    cut=tmp.weight[max.index]
    outlier.slope=((max(tmp.weight)-cut)/max(tmp.weight))/((length(tmp.weight)-max.index)/length(tmp.weight))
    if (outlier.slope>1){
      print (paste("outlier cut is ",outlier.slope,": >1, correcting outliers of weights",paste="")) 
      tmp.weight=ifelse(tmp.weight>cut,nthroot(tmp.weight,2), tmp.weight) #weight transform nthroot = root of 1/n, here we used square root
    } 
  }
  tmp.weight=sort(normalize(tmp.weight),decreasing=TRUE)
  myPermu<-list()
  print("Pre-calculating 1000 ES permutations") #corrected: precalculate all possible positive numbers less than the max number
  matchn=sapply(feature.list,function(x){length(intersect(names(tmp.weight),x))})
  matchn=sort(unique(matchn))
  k=0
  for(j in matchn){
    tmp.permu<-perm_weightedKS(tmp.weight,j,1000)
    myPermu[[j]]<-tmp.permu
    names(myPermu)[j]<-j   
    k=k+1
    if(k %% 100==0){
      print(paste("Finished permutations for ",k," matched numbers",sep=""))
    }
  }
  print("Pre-calculation finished")
  ks.result<-weightedKSV2(tmp.weight,feature.list,myPermu=myPermu,p.cut=p.cut,correct.overfit=correct.overfit,min.numOnList=minsize)
  return(ks.result)
}
#transform the negative tail of the weight into small positive values
weight.transform<-function(weight){
weight=ifelse(weight<=0,min(weight[weight>0])*((weight-min(weight[weight<=0]))/(max(weight[weight<=0])-min(weight[weight<=0]))),weight)
return(weight)
}
#run weighted KS for a weighted subject list using weightedKSV2 module

run.weightedKS<-function(weight,signed=TRUE,feature.list,minsize=10,correct.overfit=FALSE,correct.outlier=FALSE,transformNegWeight=FALSE){
  feature.select<-lapply(feature.list,filter,filter=names(weight),inlist=TRUE)
  feature.select<-feature.select[which(sapply(feature.select,function(x) length(x)>=minsize))]
  if (signed==TRUE){
    if (transformNegWeight==TRUE){
      posWeight=weight.transform(weight)
      negWeight=weight.transform(-weight)
    }else{
      posWeight=weight
      negWeight=-weight
    }
    ks.descending<-weightedKS.batch(weight=posWeight,feature.select,preCal = TRUE,high.enrich = TRUE,correct.overfit=correct.overfit,correct.outlier=correct.outlier,minsize=minsize)
    ks.ascending<-weightedKS.batch(weight=negWeight,feature.select,preCal = TRUE,high.enrich = TRUE,correct.overfit=correct.overfit,correct.outlier=correct.outlier,minsize=minsize) 
  }else{
    ks.descending<-weightedKS.batch(weight,feature.select,preCal = TRUE,high.enrich = TRUE,correct.overfit=correct.overfit,correct.outlier=correct.outlier,minsize=minsize)
    ks.ascending<-weightedKS.batch(weight,feature.select,preCal = TRUE,high.enrich = FALSE,correct.overfit=correct.overfit,correct.outlier=correct.outlier,minsize=minsize) 
  }
  ks.all=list(ks.descending=ks.descending,ks.ascending=ks.ascending)
  return(ks.all)
}

# calculate ES based on weights (sorted in decreasing order)
cal_ES<-function(weight,posList){
  if (sum(weight<0)>0) {
    stop (paste("negative values found in weight, aborting calculation of ES",weight,sep=" "))#there cannot be nagative values in weight
  }
  weight<-sort(weight,decreasing=TRUE)
  onIndice<-which(names(weight) %in% posList)
  stepDown<-1/(length(weight)-length(onIndice))
  norWeight<-rep(-stepDown,length(weight))
  names(norWeight)<-names(weight)
  norWeight[onIndice]<-weight[onIndice]/sum(weight[onIndice]) #corrected: should not be abs(weight[onIndice]/sum(weight[onIndice])) and there cannot be negative values
  cumNorWeight<-cumsum(norWeight)
  minES<-min(cumNorWeight)
  maxES<-max(cumNorWeight)
  if (abs(minES)>=abs(maxES)){
    return(0)    
  }else{
    return(maxES)    
  }
}

#calculate permutated ES
perm_weightedKS<-function(weights,numOnList,nPermu=1000){
  permu<-c()
  while(length(permu)<nPermu){
    tmp.genolist<-sample(names(weights),numOnList)
    permES<-cal_ES(weights,tmp.genolist)
    permu<-append(permu,permES)
  }
  permu<-as.numeric(permu)
  return(permu)
}



#weights,compare.list,myPermu,p.cut=0.05,correct.overfit=FALSE,min.numOnList=5
weightedKSV2<-function(weights,compare.list,myPermu,p.cut=0.05,correct.overfit=F,min.numOnList=5){ 
  cat("correct.overfit=", correct.overfit, "\n")
  mytest = lapply(compare.list, function (x) {
    ES = cal_ES(weights, x)
    onIndice = which(names(weights) %in% x)
    numOnList = length(onIndice)
    if(numOnList < min.numOnList|ES==0) {
      return(c(0, 0, "", rep("",length(weights))))
    } else {
      permu = myPermu[[numOnList]]
      NES = ES/mean(permu)
      pValue = length(permu[which(permu >= ES)])/length(permu)
      NES.subject = rep("",length(weights))
      if (pValue<p.cut & correct.overfit==T){
        for (i in onIndice){
          geneset.i=x[x != names(weights[i])]
          ESi<-cal_ES(weights,geneset.i)
          NESi<-ESi/mean(myPermu[[numOnList-1]])
          NES.subject=replace(NES.subject, i, NESi)      
        } 
      }
      return(c(NES,pValue,"NA",NES.subject))
    }
  })
  mytest.out = bind_rows(mytest) %>% t()
  colnames(mytest.out) = c("NES","pValue","qValue",names(weights))
  mytest.out = cbind(Compare.List = rownames(mytest.out), mytest.out)
  rownames(mytest.out) = NULL
  mytest.out = mytest.out[mytest.out[, 'qValue'] == "NA", ]
  mytest.out[, 'qValue']=qvalue(as.numeric(mytest.out[, 'pValue']), pi0 = 1)$qvalues
  mytest.out=mytest.out[as.numeric(mytest.out[, 'pValue'])<p.cut,]
  mytest.out = as.data.frame(mytest.out,stringsAsFactors=F)
  if (correct.overfit==T){
    return(mytest.out)
  }else{
    return(mytest.out[,1:which(colnames(mytest.out)=='qValue')])
  }
}
#perform D-CSEA2 for pathway enrichment analysis based on a dichotomous target gene list
CSEA2<-function(target.score,compare.list,p.cut=0.05,minsize=5,min.numOnList=5,transformNegWeight=FALSE){
    if (!exists("compare.list")){
      stop("please provide the list of concepts to variable: feature.list")
    }
    #tmp.weight=target.score$uniConSig
    #names(tmp.weight)=target.score$subjectID
    if (transformNegWeight==TRUE){
      target.score=weight.transform(target.score)
    }else{
      target.score=normalize(target.score)
    }
    tmp.weight=target.score
    tmp.weight=tmp.weight[which(tmp.weight!=0)]
    if (length(tmp.weight)<50){
      stop(paste("There are only",length(tmp.weight),"genes have uniConSig scores, cannot perform CSEA")) 
    }else if(length(tmp.weight)<1000){
      print (paste("There are only",length(tmp.weight),"genes have uniConSig scores, CSEA results may not be reliable"))
    }else {
      print (paste("There are",length(tmp.weight),"genes have uniConSig scores"))
    }
    compare.list=lapply(compare.list,filter,filter=names(tmp.weight),inlist=TRUE)
    compare.list<-compare.list[which(sapply(compare.list,function(x) length(x)>=minsize))]
    print("Pre-calculating 1000 ES permutations") #corrected: precalculate all possible positive numbers less than the max number
    myPermu<-list()
    matchn=sapply(compare.list,function(x){length(intersect(names(tmp.weight),x))})
    matchn=sort(unique(matchn))
    matchn=matchn[matchn>=min.numOnList]
    k=0
    for(j in matchn){
      tmp.permu<-perm_weightedKS(tmp.weight,j,2000)
      myPermu[[j]]<-tmp.permu
      names(myPermu)[j]<-j    
      k=k+1
      if(k %% 50==0){
        print(paste("Finished permutations for ",k," matched numbers",sep=""))
      }
    }
    print("Pre-calculation finished")
    ks.result<-weightedKSV2(tmp.weight,compare.list,myPermu=myPermu,p.cut=p.cut,min.numOnList=min.numOnList)
    ks.result=ks.result[order(ks.result$NES,decreasing = TRUE),]
    if (nrow(ks.result)==0){
      print("no pathway is found to be significant by KS test")
      return(NULL)
    }else{
      row.names(ks.result)=1:nrow(ks.result)
      return(ks.result)
    }
}

benchmark.uniConSig<-function(tmp.weight,training.list,testing.list=NA){
  if (suppressWarnings(is.na(testing.list))){
    compare.list=list(Training.List=training.list)
  }else{
    compare.list=list(Training.List=training.list,Testing.List=testing.list) 
  }
  matchn=sapply(compare.list,function(x){length(intersect(names(tmp.weight),x))})
  matchn=sort(unique(matchn))
  myPermu<-list()
  k=0
  for(j in matchn){
    tmp.permu<-perm_weightedKS(tmp.weight,j,1000)
    myPermu[[j]]<-tmp.permu
    names(myPermu)[j]<-j  
    k=k+1
    if(k %% 50==0){
      print(paste("Finished permutations for ",k," matched numbers",sep=""))
    }
  }
  ks.result<-weightedKSV2(tmp.weight,compare.list,myPermu=myPermu,p.cut=1)  
  return(ks.result)
}

foldCrossVal.benchmarkUniConSig<-function(target.list,feaure.list,preCalmatrix,fold=5){
  testSets=foldCrossGroups(target.list,fold=fold)
  ca.uniConSig=list()
  benchmark=data.frame(matrix(ncol = 4, nrow = 0))
  for (i in 1:length(testSets)){
    trainSet=target.list[which(!(target.list %in% testSets[[i]]))]
    ca.uniConSig[i]=list(i=cal.uniConSig(target.list=trainSet,feature.list,preCalmatrix,rm.overfit=T))
    tmp.weight=ca.uniConSig[[i]]$uniConSig
    names(tmp.weight)=ca.uniConSig[[i]]$subjectID
    benchmark=rbind(benchmark,benchmark.uniConSig(tmp.weight,training.list=trainSet,testing.list = testSets[[i]]))
  }
  benchmark$permutation=rep(1:fold,each=2)
  benchmark.tidy= benchmark %>% spread (Compare.List,NES)
  return(benchmark.tidy)
}
#disambiguate pathway result
disambiguation.CSEA<-function(GEA.result,uniConSig.result,compare.list,upPathways=NULL,topn=30,p.cut=0.01){ #upPathways TRUE or FALSE (for W-CSEA) or NULL (NULL is for D-CSEA)
  topPathway=GEA.result$Compare.List[1:topn]
  topPathway.list=compare.list[topPathway]
  topPathway.del=list()
  for (i in 1:length(topPathway.list)){
    pathwayx=topPathway.list[[i]]
    for (j in 1:length(topPathway.list)){
      pathwayy=topPathway.list[[j]]
      if (i==j){
        next
      }
      list=list(a=pathwayx[which(!pathwayx %in% pathwayy)])
      names(list)=c(paste(names(topPathway.list)[i],":",names(topPathway.list)[j],sep=""))
      topPathway.del=c(topPathway.del,list)
    }
  }
  if (is.null(upPathways)){
    CSEA.del<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$uniConSig), uniConSig.result$subjectID),compare.list=topPathway.del,p.cut=0.05) 
  }else if (upPathways==TRUE){
    CSEA.del<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$up.uniConSig), uniConSig.result$subjectID),compare.list=topPathway.del,p.cut=0.05)
  }else{
    CSEA.del<-CSEA2(target.score=setNames(as.numeric(uniConSig.result$down.uniConSig), uniConSig.result$subjectID),compare.list=topPathway.del,p.cut=0.05)
  } 
  tmp.bin <- matrix(0, nrow = length(topPathway), ncol = length(topPathway))
  dimnames(tmp.bin) <- list(topPathway, topPathway)
  ij=which(tmp.bin==0, arr.ind = T)
  for (k in 1:nrow(ij)){
    delname=paste(rownames(tmp.bin)[ij[k,"row"]],colnames(tmp.bin)[ij[k,"col"]],sep=":")
    if (delname %in% CSEA.del$Compare.List){
      tmp.bin[ij[k,"row"],ij[k,"col"]]=as.numeric(CSEA.del$pValue[CSEA.del$Compare.List==delname])
    }else{
      tmp.bin[ij[k,"row"],ij[k,"col"]]=1
    }
  }
  index=which(tmp.bin>p.cut,arr.ind = TRUE)
  remove=c()
  for (x in 1:nrow(index)){
    rowname=rownames(tmp.bin)[index[x,"row"]]
    colname=colnames(tmp.bin)[index[x,"col"]]
    if (rowname==colname | rowname %in%remove |colname %in% remove ){
      next
    }
    remove=c(remove,rowname)
  }
  filter.result=GEA.result[1:topn,]
  filter.result=filter.result[which(!filter.result$Compare.List %in% remove),]
  return(list(pathway.filtered=filter.result,pathway.todel=remove,disambiguate.matrix=tmp.bin))
}

#calculation the functional associations of the pathway enrichment result based on D-CSEA2 and return a similarity matrix between top N pathways
pathwayAssociation<-function(topPathway,compare.list,feature.list,preCalmatrix,minsize=10,rm.overfit=FALSE){ #upPathways TRUE or FALSE
  topPathway=topPathway[which(sapply(topPathway,function(x) length(compare.list[[x]])>=minsize))]
  topPathway.list=compare.list[topPathway]
  topPathway.CSEA=list()
  for (i in 1:length(topPathway.list)){
    pathwayx=topPathway.list[[i]]
    tmp.uniConSig=cal.uniConSig(target.list=pathwayx,feature.list=feature.list,preCalmatrix=preCalmatrix,minsize=minsize,rm.overfit=rm.overfit)
    tmp.CSEA=CSEA2(target.score=setNames(as.numeric(tmp.uniConSig$uniConSig), tmp.uniConSig$subjectID),compare.list=topPathway.list,p.cut=2,minsize=0) #make sure to use p.cut=2 and minsize=0 to obtain CSEA results for all pathways.
    topPathway.missed=names(topPathway.list)[!names(topPathway.list) %in% tmp.CSEA$Compare.List]
    if (length(topPathway.missed)>0){
      for (pathway in topPathway.missed){
        tmp.CSEA=rbind(tmp.CSEA,c(pathway,0,1,1))
      }
    }
    list=list(a=tmp.CSEA)
    names(list)=c(paste(names(topPathway.list)[i],sep=""))
    topPathway.CSEA=c(topPathway.CSEA,list)
  }
  tmp.bin <- matrix(0, nrow = length(topPathway), ncol = length(topPathway))
  dimnames(tmp.bin) <- list(topPathway, topPathway)
  ij=which(tmp.bin==0, arr.ind = T)
  for (k in 1:nrow(ij)){
    pathwayi=rownames(tmp.bin)[ij[k,"row"]]
    pathwayj=colnames(tmp.bin)[ij[k,"col"]]
    tmp.CSEA=topPathway.CSEA[[pathwayi]]
    if (pathwayi != pathwayj){
      tmp.bin[ij[k,"row"],ij[k,"col"]]=tmp.CSEA$NES[tmp.CSEA$Compare.List==pathwayj] 
    }
    #print (k)
  }
  mode(tmp.bin)<-"numeric"
  diag(tmp.bin) <- max(tmp.bin)
  return(tmp.bin)
}

pathway.corplot<-function(corMatrix,fontsize=0.8){
  corMatrix=2*(corMatrix-min(corMatrix))/(max(corMatrix)-min(corMatrix))-1
  rownames(corMatrix)=gsub("_"," ",rownames(corMatrix))
  colnames(corMatrix)=gsub("_"," ",colnames(corMatrix))
  mycolor=colorRampPalette(c("dodger blue", "white", "red"))(n = 100)
  corrplot(corMatrix, method = "circle",type ="lower",tl.col="black",col=mycolor,tl.cex=fontsize,tl.srt=70)
}
library(pheatmap)
pathway.heatmap<-function(matrixData,clustering=TRUE,distanceMatric="euclidean",fontSize=0.7,sepWidth=0.1){
  rownames(matrixData)=gsub("_"," ",rownames(matrixData))
  colnames(matrixData)=gsub("_"," ",colnames(matrixData))
  my_palette <- colorRampPalette(c("dodger blue", "white", "red"))(n = 1000)
  pheatmap(matrixData, cluster_rows=clustering, cluster_cols=clustering,color=my_palette, show_rownames = T, show_colnames = T,fontsize=fontSize,legend=T) 
  return(pathway.heatmap)
}

gene.heatmap<-function(matrixData,fontSize=0.7,sepWidth=0.1){
    my_palette <- colorRampPalette(c("dodger blue", "white", "red"))(n = 1000)
    heatmap.2(matrixData, distfun=function(matrixData) dist(matrixData, method=distanceMatric),
              col=my_palette,       # use on color palette defined earlier
              symm=F,symkey=F,symbreaks=F, 
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              Colv=FALSE,  #whether or not perform column clustering (ordering)
              Rowv=FALSE,  #whether or not perform row clustering (ordering)
              hclustfun=hclust,
              dendrogram="none",  # "row", "both"
              cexRow=fontSize, cexCol=fontSize,keysize=1,
              offsetRow = 0,
              offsetCol = 0,
              margins = c(20, 25)
              #sepwidth=c(sepWidth,sepWidth),
              #sepcolor="grey",
              #colsep=1:ncol(matrixData),
              #rowsep=1:nrow(matrixData)
    )
}

gene.heatmap2<-function(genomatrix,fontsize){
  library(matrixStats)
  library(pheatmap)
  quantile.range <- quantile(genomatrix, probs = seq(0, 1, 0.01))
  myBreaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
  myColor  <- colorRampPalette(c("skyblue", "white", "red"))(length(myBreaks) - 1)
  pheatmap(genomatrix, cluster_rows=FALSE, cluster_cols=FALSE,color=myColor, breaks=myBreaks,show_rownames = T, show_colnames = F,fontsize_row=fontsize,legend=F) 
}

draw.pathway <- function(weight,pathway,pathway.list){
  library(pacman)
  pacman::p_load(readxl, ggplot2,ggpubr,magrittr)
  
  weight<-sort(weight,decreasing=TRUE)
  figData=as.data.frame(cbind(weight,gene=names(weight)))
  figData$inPathway=ifelse(figData$gene %in% pathway.list[[pathway]],1,0)
  onIndice<-which(names(weight) %in% pathway.list[[pathway]])
  stepDown<-1/(length(weight)-length(onIndice))
  stepUp<-weight[onIndice]/sum(weight[onIndice])
  figData$norWeight<-ifelse(figData$inPathway==1,stepUp,-stepDown)
  figData$cumNorWeight<-cumsum(figData$norWeight)
  #enrichFig<-plot(figData$cumNorWeight,type="l", xlab="",ylab="" )
  
  ######### by Sanghoon##########################
  dataset<-data.frame(V1=rownames(figData),V2=figData[,1],stringsAsFactors=F)
  colnames(dataset)[1]<-"NAME"
  dataset <- dataset[match(unique(dataset$NAME),dataset$NAME),]
  dataset = subset(dataset, select=-c(NAME))
  gene.labels <- row.names(dataset)
  A <- data.matrix(dataset)
  cols <- length(A[1, ])
  rows <- length(A[, 1])
  
  # Read input gene set database
  temp<-pathway.list
  max.Ng <- length(temp)
  temp.size.G <- vector(length = max.Ng, mode = "numeric")
  for (i in 1:max.Ng) {
    temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\\t"))) - 2
  }
  
  max.size.G <- max(temp.size.G)
  min.size.G <- min(temp.size.G)
  gs <- matrix(rep("null", max.Ng * max.size.G), nrow = max.Ng, ncol = max.size.G)
  temp.names <- vector(length = max.Ng, mode = "character")
  temp.desc <- vector(length = max.Ng, mode = "character")
  gs.count <- 1
  for (i in 1:max.Ng) {
    gene.set.size <- length(unlist(strsplit(temp[[i]], "\\t")))-2
    gs.line <- noquote(unlist(strsplit(temp[[i]], "\\t")))
    gene.set.name<-gs.line[1]
    gene.set.desc<-gs.line[2]
    gene.set.tags<-vector(length = gene.set.size, mode = "character")
    for (j in 1:gene.set.size) {
      gene.set.tags[j]<-gs.line[j + 2]
    }
    existing.set <- is.element(gene.set.tags, gene.labels)
    set.size <- length(existing.set[existing.set == T])
    temp.size.G[gs.count]<-set.size
    gs[gs.count, ]<-c(gene.set.tags[existing.set], rep("null", max.size.G-temp.size.G[gs.count]))
    temp.names[gs.count]<-gene.set.name
    temp.desc[gs.count]<-gene.set.desc
    gs.count <- gs.count + 1
  }
  
  Ng <- gs.count - 1
  gs.names <- vector(length = Ng, mode = "character")
  gs.desc <- vector(length = Ng, mode = "character")
  size.G <- vector(length = Ng, mode = "numeric")
  gs.names <- temp.names[1:Ng]
  gs.desc <- temp.desc[1:Ng]
  size.G <- temp.size.G[1:Ng]
  
  N <- length(A[, 1])  # Number of Genes in weight
  Ns <- length(A[1, ])
  
  print(c("Number of genes:", N))
  print(c("Number of Gene Sets:", Ng))
  print(c("Number of samples:", Ns))
  print(c("Original number of Gene Sets:", max.Ng))
  print(c("Maximum gene set size:", max.size.G))
  print(c("Minimum gene set size:", min.size.G))
  
  Obs.indicator <- matrix(nrow = Ng, ncol = N)
  Obs.indicator<-figData$inPathway
  Obs.RES <- matrix(nrow = Ng, ncol = N)
  
  min.RES<-min(figData$cumNorWeight)
  max.RES<-max(figData$cumNorWeight)
  if (max.RES < 0.3) 
    max.RES <- 0.3
  if (min.RES > -0.3) 
    min.RES <- -0.3
  delta<-(max.RES-min.RES)*0.5
  min.plot<-min.RES-2*delta
  max.plot<-max.RES
  
  obs.s2n.fac<-unname(figData$weight)
  obs.s2n<-as.numeric(levels(obs.s2n.fac))[obs.s2n.fac]
  max.corr<-max(obs.s2n)
  min.corr<-min(obs.s2n)
  Obs.correl.vector.norm<-(obs.s2n-min.corr)/(max.corr-min.corr)*1.25*delta+min.plot
  zero.corr.line<-(-min.corr/(max.corr-min.corr))*1.25*delta+min.plot
  
  sub.string<-paste("Number of genes: n=",N, " (",length(onIndice), " in ",pathway,")",sep="", collapse="")
  main.string<-paste("Gene Set ",":", pathway)
  
  plot(1:N, figData$cumNorWeight, sub=sub.string,xlab="Rank in Ordered Dataset",
       ylab = "Running Enrichment Score (RES)",type="l",xlim = c(1,N), ylim=c(min.plot,max.plot),lwd=2,cex=1,col=2)
  
  # Ranked list metric (Signal2Noise)
  for (j in seq(1, N, 20)) {
    lines(c(j, j), c(zero.corr.line,Obs.correl.vector.norm[j]), lwd=1,cex=1,col=colors()[12])  # shading of correlation plot
  }
  
  lines(c(1, N),c(0, 0),lwd=1,lty=2,cex=1,col=1)  # zero RES line
  lines(c(which.max(figData$cumNorWeight),which.max(figData$cumNorWeight)),c(min.plot,max.plot),lwd=1,lty=3,cex=1,col=2)    # max enrichment vertical line
  
  for (j in 1:N) {
    if (Obs.indicator[j]==1) {
      lines(c(j,j), c(min.plot+1.25*delta,min.plot+1.75*delta),lwd=1,lty=1,cex=1,col=4)  # enrichment tags
    }
  }
  
  lines(1:N, Obs.correl.vector.norm, type="l",lwd=1,cex=1,col=1)
  lines(c(1,N), c(zero.corr.line,zero.corr.line),lwd=1,lty=1,cex=1,col=1)  # zero correlation horizontal line
  
  temp <- order(abs(obs.s2n), decreasing=T)   
  arg.correl<-temp[N]  ## gene position of minimum "Rank Metric Score" 
  lines(c(arg.correl,arg.correl), c(min.plot,max.plot), lwd=1,lty=3,cex=1,col=3)  # zero crossing correlation vertical line
  
  adjx <- ifelse(max.RES>0, 0, 1)
  leg.txt <- paste("Peak at ", which.max(figData$cumNorWeight), sep = "", collapse = "")
  text(x=which.max(figData$cumNorWeight), y=min.plot+1.8*delta, adj = c(adjx,0), labels = leg.txt, cex = 1)
  
  ## The gene crossing zero Rank Metric Score (the gene of the minimum of weight)
  leg.txt <- paste("Zero crossing at ", arg.correl, sep = "", collapse = "")
  text(x = arg.correl, y = min.plot + 1.95 * delta, adj = c(adjx,0), labels = leg.txt, cex = 1)
} # end of draw.pathway function

#merge upregulated and downregulated pathways. For pathways that are both up or down regulated, the most significant rank will be retained
merge.pathway<-function(up.pathway,down.pathway){
  pathway.ambiguous=intersect(up.pathway$Compare.List,down.pathway$Compare.List)
  up.pathway$Type="UP"
  down.pathway$Type="DOWN"
  if (length(pathway.ambiguous)==0){
    pathway.combine=rbind(up.pathway,down.pathway)
  }else{
    ambiguous.uprank=setNames(rownames(up.pathway)[match(pathway.ambiguous,up.pathway$Compare.List)],pathway.ambiguous)
    ambiguous.downrank=setNames(rownames(down.pathway)[match(pathway.ambiguous,down.pathway$Compare.List)],pathway.ambiguous)
    pathway.ambiguous.type=setNames(ifelse(as.numeric(ambiguous.uprank[pathway.ambiguous])<=as.numeric(ambiguous.downrank[pathway.ambiguous]),"UP","DOWN"),pathway.ambiguous)
    up.pathway=up.pathway[!up.pathway$Compare.List %in% names(pathway.ambiguous.type)[pathway.ambiguous.type=="DOWN"],]
    down.pathway=down.pathway[!down.pathway$Compare.List %in% names(pathway.ambiguous.type)[pathway.ambiguous.type=="UP"],]
    pathway.combine=rbind(up.pathway,down.pathway)
  }
  return(pathway.combine)
}
#draw network figure for significant pathways based on association file
draw.network<-function(pathway.out,assoc,NES.cut=2,node.size=2,line.thickness=1.5){
  links=matrixConvert(assoc, colname = c("from", "to", "similarity"))
  links=links[links$similarity>NES.cut,]
  if ("Type" %in% colnames(pathway.out)){
    nodes=data.frame(id=row.names(assoc),size=node.size*as.numeric(pathway.out$NES[match(row.names(assoc),pathway.out$Compare.List)]),type=pathway.out$Type[match(row.names(assoc),pathway.out$Compare.List)],stringsAsFactors = F)
  }else{
    nodes=data.frame(id=row.names(assoc),size=node.size*as.numeric(pathway.out$NES[match(row.names(assoc),pathway.out$Compare.List)]),type="NA",stringsAsFactors = F)
  }
  #https://kateto.net/netscix2016.html
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=F)  
  colrs <- c("darkolivegreen2","peachpuff")
  V(net)$color <- colrs[as.factor(V(net)$type)]#set node color
  V(net)$size <- V(net)$size #set node size
  V(net)$label <- gsub("HALLMARK_|REACTOME_|KEGG_|WP_|BIOCARTA_","",names(V(net)))
  V(net)$label <- gsub("_"," ",V(net)$label)
  E(net)$width <- E(net)$similarity*line.thickness
  E(net)$arrow.size <-0
  E(net)$color <- 'cornsilk2'
  #dev.new(width=12, height=5, unit="in")
  windowsFonts("Arial" = windowsFont("Arial"))
  plot(net,vertex.label.cex=.7,vertex.label.color="black",vertex.label.cex = 0.4, vertex.label.degree = 4,vertex.frame.color="ivory4",asp=0.5,layout=layout_with_fr)
}


