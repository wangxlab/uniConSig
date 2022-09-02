library(stringr)
library(vegan)
library(qvalue)
library(tidyr)
library(moments)
library(limma)
library(dplyr)
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

#read genotype gmt file into a list; genotype name and class are seperated with "#"
read_gmt<-function(fileName,min=5){
  con=file(fileName,"r")
  i<-0
  tmp.posList<-list()
  while (length(line<-readLines(con,n=1))>0){
    if (grepl("^#", line)){
      #print(paste(line,"skipped",sep=" "))
      next
    }
    tmp.line <- unlist(strsplit(line,"\t"))
    tmp.name <- tmp.line[1]
    if (length(tmp.line)>=min+2){
      i=i+1
      tmp.posList[[i]] <- tmp.line[3:length(tmp.line)]
      names(tmp.posList)[i]<-tmp.name  
      if(i %% 2000==0){
        print(paste("Loaded ",i," features",sep=""))
      }
    }
  }
  close(con)
  return(tmp.posList)
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

#precalculate feature redundancy
batch_calSimilarity<-function(feature.list,feature.preCalfile,cutoff=0.1,minsize=20,method){
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

#calculate uniConSig scores
cal.uniConSig<-function(target.list,feature.list,preCalmatrix,minsize=20,weight.cut=0.05,power=1,root=0.5,ECNpenalty=1,method="Ochiai",rm.overfit=TRUE){ #method=="Ochiai" or "Jaccard"
  target.weight=batch_CalWeight(target.list,feature.list,method=method)
  target.weight=target.weight[as.numeric(target.weight$Compare.Size)>minsize&as.numeric(target.weight$Weight)>weight.cut,]
  if (nrow(target.weight)<10){
    stop(paste("There are only",nrow(target.weight),"signature features with weights, this means that the target list is functionally heterogeneous, please lower the cutoff for weight.cut"))
  }
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
    ECN=sum(1/(tmp.data$Epsilon^root))
    uniConSig<-sum((tmp.data$Weight4Cal^power)/(tmp.data$Epsilon^root))/(ECN^ECNpenalty)
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
#calculated weighted KS using weightedKSV2 based on a gene list with t-values as weights (such as the gene list sorted by t-values)
weightedKS.batch<-function(weight,feature.list,preCal=TRUE,high.enrich=TRUE,p.cut=0.05,correct.overfit=FALSE,correct.outlier=FALSE,minsize=10){
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
  ks.result<-weightedKSV2(tmp.weight,feature.list,myPermu=myPermu,p.cut=p.cut,correct.overfit=correct.overfit)
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
    ks.descending<-weightedKS.batch(weight=posWeight,feature.select,preCal = TRUE,high.enrich = TRUE,correct.overfit=correct.overfit,correct.outlier=correct.outlier)
    ks.ascending<-weightedKS.batch(weight=negWeight,feature.select,preCal = TRUE,high.enrich = TRUE,correct.overfit=correct.overfit,correct.outlier=correct.outlier) 
  }else{
    ks.descending<-weightedKS.batch(weight,feature.select,preCal = TRUE,high.enrich = TRUE,correct.overfit=correct.overfit,correct.outlier=correct.outlier)
    ks.ascending<-weightedKS.batch(weight,feature.select,preCal = TRUE,high.enrich = FALSE,correct.overfit=correct.overfit,correct.outlier=correct.outlier) 
  }
  ks.all=list(ks.descending=ks.descending,ks.ascending=ks.ascending)
  return(ks.all)
}
####Calculate uniConSig scores based on ks test results
cal.uniConSig.ks<-function(up.ks,down.ks,preCalmatrix,feature.list,outfile,p.cut=0.01,q.cut=0.25,NES.cut=0,power=1,root=1,ECNpenalty=0.5,correct.overfit=FALSE){
  up.ks[is.na(up.ks)] <- 0
  down.ks[is.na(down.ks)] <- 0
  up.ks.sub<-up.ks[up.ks$NES>NES.cut & up.ks$pValue<p.cut & up.ks$qValue<q.cut,-c(3,4)]
  up.ks.mat <- as.matrix(up.ks.sub[-1])
  row.names(up.ks.mat) <- up.ks.sub$Compare.List
  down.ks.sub<-down.ks[down.ks$NES>NES.cut & down.ks$pValue<p.cut & down.ks$qValue<q.cut,-c(3,4)]
  down.ks.mat <- as.matrix(down.ks.sub[-1])
  row.names(down.ks.mat) <- down.ks.sub$Compare.List
  up.list=list()
  down.list=list()
  result<-data.frame(matrix(ncol = 3, nrow = 0))
  for (i in 1:nrow(preCalmatrix)){
    tmp.line<-unlist(strsplit(preCalmatrix[i,],"\t"))
    subjectID=as.character(tmp.line[1])
    if(length(tmp.line)==2){
      next
    }
    tmp.epsilon<-as.data.frame(do.call(rbind, strsplit(tmp.line[3:length(tmp.line)],"@")))
    colnames(tmp.epsilon)<-c("Compare.List","Epsilon")
    ECN=sum(1/(as.numeric(tmp.epsilon[,"Epsilon"])^root))
    up.epsilon<-tmp.epsilon[tmp.epsilon$Compare.List %in% row.names(up.ks.mat),]
    down.epsilon<-tmp.epsilon[tmp.epsilon$Compare.List %in% row.names(down.ks.mat),]
    if (correct.overfit==TRUE){
      if (subjectID %in% names(up.ks)){
        up.weight<- up.ks.mat[which(row.names(up.ks.mat) %in% as.character(up.epsilon$Compare.List)),c("NES",subjectID),drop=F]
        up.weight <- as.data.frame(replace(up.weight[,"NES"], which(up.weight[,subjectID]!=0),up.weight[which(up.weight[,subjectID]!=0),subjectID]))
      }else{
        up.weight<- as.data.frame(up.ks.mat[which(row.names(up.ks.mat) %in% up.epsilon$Compare.List),"NES"])
      }      
    }else{
      up.weight<- as.data.frame(up.ks.mat[which(row.names(up.ks.mat) %in% up.epsilon$Compare.List),"NES"])
    }
    up.weight$Compare.List <- as.character(row.names(up.weight))
    colnames(up.weight)=c("NES","Compare.List")
    up.data<-as.matrix(merge(up.epsilon,up.weight,by.x="Compare.List",by.y="Compare.List"))
    up.data=cbind(up.data,normWeight=(as.numeric(up.data[,"NES"])^power)/(as.numeric(up.data[,"Epsilon"])^root))
    score.sensitive<-sum(as.numeric(up.data[,"normWeight"]))/(ECN^ECNpenalty)
    up.data=as.data.frame(up.data)
    colnames(up.data)[-1]=c(paste(subjectID,colnames(up.data[-1]),sep=":"))
    up.list[[paste(subjectID,sep="")]]<-up.data
    if (correct.overfit==TRUE){    
      if (paste("X",subjectID,sep="") %in% names(down.ks)){
        down.weight<- down.ks.mat[which(row.names(down.ks.mat) %in% down.epsilon$Compare.List),c("NES",paste("X",subjectID,sep="")),drop=F]
        down.weight <- as.data.frame(replace(down.weight[,"NES"], which(down.weight[,paste("X",subjectID,sep="")]!=0),down.weight[which(down.weight[,paste("X",subjectID,sep="")]!=0),paste("X",subjectID,sep="")]))
      }else{
        down.weight<- as.data.frame(down.ks.mat[which(row.names(down.ks.mat) %in% down.epsilon$Compare.List),"NES"])
      }
    }else{
      down.weight<- as.data.frame(down.ks.mat[which(row.names(down.ks.mat) %in% down.epsilon$Compare.List),"NES"])
    }
    down.weight$Compare.List <- as.character(row.names(down.weight))
    colnames(down.weight)=c("NES","Compare.List")
    res.data<-as.matrix(merge(down.epsilon,down.weight,by.x="Compare.List",by.y="Compare.List"))
    res.data=cbind(res.data,normWeight=(as.numeric(res.data[,"NES"])^power)/(as.numeric(res.data[,"Epsilon"])^root))
    score.resistant<-sum(as.numeric(res.data[,"normWeight"]))/(ECN^ECNpenalty)
    res.data=as.data.frame(res.data)
    colnames(res.data)[-1]=c(paste(subjectID,colnames(res.data[-1]),sep=":"))
    down.list[[paste(subjectID,sep="")]]<-res.data
    result[nrow(result) + 1,]<-c(tmp.line[1],score.sensitive,score.resistant)
    if(i %% 2000==0){
      print(paste("Processed ",i," subjects",sep=""))
    }
  }
  colnames(result)<-c("subjectID","up.uniConSig","down.uniConSig")
  return(result)
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


#calculated weighted KS and output NES, p-value, q-value 
weightedKSV2<-function(weights,compare.list,myPermu,p.cut=0.05,correct.overfit=FALSE,min.numOnList=5){ 
  KS.out<-data.frame(matrix(ncol = 4, nrow = 0))
  for(k in 1:length(compare.list)){
    if(k %% 1000==0){
      print(paste("Finished ",k," features",sep=""))
    }
    ES<-cal_ES(weights,compare.list[[k]])
    onIndice<-which(names(weights) %in% compare.list[[k]])
    numOnList=length(onIndice)
    if(numOnList<min.numOnList|ES==0){ #required at least 5 index cell lines for each genotype or set NES to zero
      KS.out[nrow(KS.out) + 1,]<-c(names(compare.list[k]),0,1,"NA")
    }else{
      permu<-myPermu[[numOnList]]
      NES<-ES/mean(permu)
      names(NES)<-names(compare.list[k])
      pValue<-length(permu[which(permu>=ES)])/length(permu)
      names(pValue)<-names(compare.list[k])
      KS.out[nrow(KS.out) + 1,]<-c(names(compare.list[k]),NES,pValue,"NA")
    }
  }
  colnames(KS.out)=c("Compare.List","NES","pValue","qValue")
  if (length(KS.out$pValue)<2){
    KS.out$qValue=rep("NA",nrow(KS.out))
  }else{
    KS.out$qValue=qvalue(as.numeric(KS.out$pValue), pi0 = 1)$qvalues 
  }
  KS.out=KS.out[KS.out['pValue']<p.cut,]
  if (correct.overfit==TRUE){
    print ("Correcting NES of index subjects")
    result<-data.frame(matrix(ncol = 4 + length(weights), nrow = 0))
    for (k in 1:nrow(KS.out)){
      feature.k=feature.list[[KS.out[k,"Compare.List"]]]
      onIndice=which(names(weights) %in% feature.k)
      numOnList=length(onIndice)
      NES.subject=rep("",length(weights))
      for (i in onIndice){
        feature.i=feature.k[feature.k != names(weights[i])]
        ESi<-cal_ES(weights,feature.i)
        NESi<-ESi/mean(myPermu[[numOnList-1]])
        names(NESi)=names(weights[i])
        NES.subject=replace(NES.subject, i, NESi)      
      }
      result[nrow(result) + 1,]<-c(KS.out[k,],NES.subject)
    }
    colnames(result)=c("Compare.List","NES","pValue","qValue",names(weights))
    return(result)
  }else{
    return(KS.out)
  }
}
#perform D-CSEA2 for pathway enrichment analysis based on a dichotomous target gene list
CSEA2<-function(target.score,compare.list,p.cut=0.05,minsize=5,min.numOnList=5,transformNegWeight=FALSE){
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
      tmp.permu<-perm_weightedKS(tmp.weight,j,1000)
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
    row.names(ks.result)=1:nrow(ks.result)
    return(ks.result)
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

foldCrossVal.benchmarkUniConSig<-function(target.list,feaure.list,preCalmatrix,fold=5,minsize=5,weight.cut=0,power=1,root=1,ECNpenalty=0.5,method="Ochiai"){
  testSets=foldCrossGroups(target.list,fold=fold)
  ca.uniConSig=list()
  benchmark=data.frame(matrix(ncol = 4, nrow = 0))
  for (i in 1:length(testSets)){
    trainSet=target.list[which(!(target.list %in% testSets[[i]]))]
    ca.uniConSig[i]=list(i=cal.uniConSig(target.list=trainSet,feature.list,preCalmatrix,minsize=minsize,weight.cut=weight.cut,power=power,root=root,ECNpenalty=ECNpenalty,method=method))
    tmp.weight=ca.uniConSig[[i]]$uniConSig
    names(tmp.weight)=ca.uniConSig[[i]]$subjectID
    benchmark=rbind(benchmark,benchmark.uniConSig(tmp.weight,training.list=trainSet,testing.list = testSets[[i]]))
  }
  benchmark$permutation=rep(1:fold,each=2)
  benchmark.tidy= benchmark %>% spread (Compare.List,NES)
  return(benchmark.tidy)
}
#disambiguate pathway result
disambiguation.CSEA<-function(GEA.result,uniConSig.result,compare.list,upPathways,topn=30,p.cut=0.01){ #upPathways TRUE or FALSE (for W-CSEA) or NULL (NULL is for D-CSEA)
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
disambiguation.GSEA<-function(GEA.result,weight,compare.list,topn=30,p.cut=0.01,transformNegWeight=FALSE){ #upPathways TRUE or FALSE (for W-GSEA) or NULL (NULL is for D-GSEA)
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
  GSEA.del<-CSEA2(target.score=weight,compare.list=topPathway.del,p.cut=0.05,minsize=5,transformNegWeight=transformNegWeight)
  tmp.bin <- matrix(0, nrow = length(topPathway), ncol = length(topPathway))
  dimnames(tmp.bin) <- list(topPathway, topPathway)
  ij=which(tmp.bin==0, arr.ind = T)
  for (k in 1:nrow(ij)){
    delname=paste(rownames(tmp.bin)[ij[k,"row"]],colnames(tmp.bin)[ij[k,"col"]],sep=":")
    if (delname %in% GSEA.del$Compare.List){
      tmp.bin[ij[k,"row"],ij[k,"col"]]=as.numeric(GSEA.del$pValue[GSEA.del$Compare.List==delname])
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
pathwayAssociation<-function(topPathway,compare.list,feature.list,preCalmatrix,topn=3,minsize=10,rm.overfit=FALSE){ #upPathways TRUE or FALSE
  topPathway=topPathway[which(sapply(topPathway,function(x) length(compare.list[[x]])>=minsize))]
  topPathway.list=compare.list[topPathway]
  topPathway.CSEA=list()
  for (i in 1:length(topPathway.list)){
    pathwayx=topPathway.list[[i]]
    tmp.uniConSig=cal.uniConSig(target.list=pathwayx,feature.list=feature.list,preCalmatrix=preCalmatrix,minsize=minsize,weight.cut=0.05,power=1,root=1,ECNpenalty=0.5,method="Ochiai",rm.overfit=rm.overfit)
    tmp.CSEA=CSEA2(target.score=setNames(as.numeric(tmp.uniConSig$uniConSig), tmp.uniConSig$subjectID),compare.list=topPathway.list,p.cut=2,minsize=0) #make sure to use p.cut=2 and minsize=0 to obtain CSEA results for all pathways.
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
    print (k)
  }
  mode(tmp.bin)<-"numeric"
  diag(tmp.bin) <- max(tmp.bin)
  return(tmp.bin)
}
#perform Limma analysis of gene expression data comparing two groups
limmaDGE<-function(gctFile,clsFile){ # the comparison is based on increasingly sorted orders of group labels with positive t value showing upregulated in class 1
  expData=read.table(gctFile, stringsAsFactors=F, skip=2, header=T, row.names=NULL, check.names=F, fill=TRUE,sep="\t")
  expUnique<-expData[which(!duplicated(expData[,1])),]
  rownames(expUnique)<-expUnique[,1]
  expUnique=expUnique[,-c(1,2)]
  expUnique=data.frame(lapply(expUnique, function(x) as.numeric(as.character(x))),
                        check.names=F, row.names = rownames(expUnique))
  expUnique<-expUnique[(which(!apply(expUnique,1,var)==0)),]  ## Remove genes (rows) that have zero variance
  myClass=read.table(clsFile, sep="\t", stringsAsFactors=F, header=F, skip=1, row.names=NULL, check.names=F, fill=TRUE)
  myClass<-as.factor(myClass[1,])
  label=sort(unique(as.character(myClass)))
  fit<-lmFit(expUnique,design=model.matrix(~myClass))
  fit.eBayes <- eBayes(fit)
  limmaOut<-topTable(fit.eBayes,number=100000000, adjust.method="BH")  # BH, Benjamini and Hochberg FDR adjusted p-value.
  colnames(limmaOut)<-c("Log2FC","AveExpr", "T.Value","P.Value","Q.Value","DGE.odds")
  limmaOut$Signed.Q.Value=-log10(limmaOut$Q.Value)*sign(limmaOut$T.Value)
  limmaOut<-limmaOut[order(limmaOut$Signed.Q.Value,decreasing = TRUE),]
  return(limmaOut)
}
library(gplots)
library(RColorBrewer)
pathway.heatmap<-function(matrixData,clustering=TRUE,distanceMatric="euclidean",fontSize=0.7,sepWidth=0.1){
  my_palette <- colorRampPalette(c("dodger blue", "white", "red"))(n = 1000)
  if (clustering==TRUE){
    heatmap.2(matrixData, distfun=function(matrixData) dist(matrixData, method=distanceMatric),
                       col=my_palette,       # use on color palette defined earlier
                       symm=F,symkey=F,symbreaks=F, 
                       density.info="none",  # turns off density plot inside color legend
                       trace="none",         # turns off trace lines inside the heat map
                       hclustfun=hclust,
                       dendrogram="none",  # "row", "both"
                       cexRow=fontSize, cexCol=fontSize,keysize=1,
                       offsetRow = 0,
                       offsetCol = 0,
                       margins = c(20, 25),
                       sepwidth=c(sepWidth,sepWidth),
                       sepcolor="grey",
                       colsep=1:ncol(matrixData),
                       rowsep=1:nrow(matrixData)
    )
  }else{
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
                       margins = c(20, 25),
                       sepwidth=c(sepWidth,sepWidth),
                       sepcolor="grey",
                       colsep=1:ncol(matrixData),
                       rowsep=1:nrow(matrixData)
    )
  }
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


draw.pathway <- function(weight,pathway,pathway.list){
  library(readxl)
  library(ggplot2)
  library(ggpubr)
  weight<-sort(weight,decreasing=TRUE)
  figData=as.data.frame(cbind(weight,gene=names(weight)))
  figData$inPathway=ifelse(figData$gene %in% pathway.list[[pathway]],1,0)
  onIndice<-which(names(weight) %in% pathway.list[[pathway]])
  stepDown<-1/(length(weight)-length(onIndice))
  stepUp<-weight[onIndice]/sum(weight[onIndice])
  figData$norWeight<-ifelse(figData$inPathway==1,stepUp,-stepDown)
  figData$cumNorWeight<-cumsum(figData$norWeight)
  enrichFig<-plot(figData$cumNorWeight,type="l")
  return(enrichFig)
}

draw.pathway.v2 <- function(weight,pathway,pathway.list){
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
  
  ######### I added from this line##########################
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
  
  sub.string<-paste("Number of genes: ",N, " (in list), ",length(onIndice), " (in gene set)",sep="", collapse="")
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

zlim = function(x, lim) { 
  x[which(x < min(lim))] = min(lim)
  x[which(x > max(lim))] = max(lim)
  return(x)
}

## Make columnColor function to map metadata category to color
columnColor<-function(annotations){
  colorsVector=ifelse(annotations["StudyType"]==label[1], "dodgerblue1", 
                      ifelse(annotations["StudyType"]==label[2], "coral", "blue"))
  return(colorsVector)
}

draw.heatmap <-function(gctFile, clsFile, weight, pathway,pathway.list) {
    library(pacman)
    pacman::p_load(gplots, RColorBrewer)
  
    myTableDF=read.table(gctFile, sep="\t", stringsAsFactors=F, header=F, row.names=NULL, check.names=F, fill=TRUE)[-1,-2]
    myTableNoDup<-myTableDF[which(!duplicated(myTableDF[,1])),]#REVISE
    colnames(myTableNoDup)<-myTableNoDup[1,]
    rownames(myTableNoDup)<-myTableNoDup[,1]
    
    myTableNumericDF<-as.data.frame(sapply(myTableNoDup[-1,-1], as.numeric))  ## convert character to numeric in data.frame
    rownames(myTableNumericDF)<-rownames(myTableNoDup[-1,-1])
    
    myTableNumericDF1<-log(myTableNumericDF+1,2)  ## getting log2
    zScoreDF<-scale(t(myTableNumericDF1),center=TRUE,scale=TRUE)  # z-score
    zScoreTransposed<-t(zScoreDF)
    myTableMtrx<-data.matrix(zScoreTransposed[(which(!apply(zScoreTransposed,1,var)==0)),])  # Remove genes (rows) that have zero variance,
    
    geneSet<-names(weight)[which(names(weight)%in% pathway.list[[pathway]])]     ### Make the expression matrix subset with geneSet
    mat_data<-myTableMtrx[which(rownames(myTableMtrx)%in% geneSet),]
    
    ## Get the phenotype class factors.
    myClassFactorDF=read.table(clsFile, sep="\t", stringsAsFactors=F, header=F, row.names=NULL, check.names=F, fill=TRUE)[-1,]
    fac<-as.factor(myClassFactorDF[1,])
    label=sort(unique(as.character(fac)))
    
    annotationData<-cbind(colnames(myTableMtrx), t(myClassFactorDF[1,]))    ## make annotation data
    annotData<-as.data.frame(annotationData)
    colnames(annotData)<-c("sampleID","StudyType")
    
    my_palette<-colorRampPalette(c("sky blue", "white","red"))(n=299)  # creates a own color palette from red to green
    sampleColors = columnColor(annotData)
    heatmap.2(zlim(mat_data, c(-2,2)),  ColSideColors=sampleColors,  ## "ColSideColors" is to make color annotation for samples
              distfun=function(mat_data) dist(mat_data, method=distanceMatric),
              labRow=NULL, cexRow=0.7, labCol=FALSE,  ## NULL will print the labels. FALSE doesn't print. 
              breaks=seq(-2,2,length.out=300),  # just 1 larger than my_palette(n=299)
              symm=F,symkey=F,symbreaks=T, # scale="none",symbreaks should be TRUE to display a big(long row) heatmap
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins =c(12,12),     # widens margins around plot   this is for pdf
              col=my_palette,       # use on color palette defined earlier
              reorderfun=function(d, w) reorder(d, w, hierarchicalTest.FUN = mean),
              Colv=FALSE,  # to not perform column clustering (ordering)
              Rowv=FALSE,  # to not order rows.
              #key=FALSE,
              dendrogram="none",  ## "row", "both", "none"
              #scale="row"  # To make color difference by samples, not by genes  #If i use "breaks" parameter, I can't use scale parameter. 
    )
}    
    
    
  
