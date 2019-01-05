###########################################################uniConSig#######################################################
#' @title calculate concept weight
#' @param conceptVec A vector containing only the element ids of a concept
#' @param trList A vector containing all the training gene ids.
#' @description analyze two gene sets(one training, one concept), return the number of intersected and united genes
#' @return "numberOfIntersection_numberOfUnion"
#' @examples
#' tmp.conceptVec<-1:10
#' trList.call<-4:12
#' tmp.data<-cal_conceptWeight(tmp.conceptVec,trList.call) #split the number of intersection/union.
#' @export
cal_conceptWeight<-function(conceptVec,trList){ #conceptVec contains only the element ids of the concept
    geneID.all<-length(union(trList,conceptVec))
    geneID.intersect<-length(intersect(trList,conceptVec))
    return(paste(geneID.intersect,"_",geneID.all,sep=""))
}
#' @title calculate concept weight
#' @param conceptVec A vector containing only the element ids of a concept
#' @param trList A vector containing all the training gene ids.
#' @description analyze two gene sets(one training, one concept), return the number of intersected and united genes
#' @return "numberOfIntersection_numberOfUnion"
#' @examples
#' tmp.conceptVec<-1:10
#' trList.call<-4:12
#' tmp.data<-cal_conceptWeight(tmp.conceptVec,trList.call) #split the number of intersection/union.
#' @export
library(vegan)
cal_conceptWeight_continuous<-function(conceptVec,triVec){ #conceptVec contains only the element ids of the concept
  geneID.all<-length(union(triVec,conceptVec))
  geneID.intersect<-length(intersect(triVec,conceptVec))
  return(paste(geneID.intersect,"_",geneID.all,sep=""))
}
#' @title load concept database from gmt, calculate concept weights, process the file one line at a time
#' @param trList A vector containing all the training gene ids.
#' @description calculate concept weight between the training gene set and concepts, with/without penalty for over-representing gene
#' @return A matrix containing two types of concept weights
#' @examples
#' trList.call<-1:100
#' concept.info<-cal_weightMatrix(trList.call)
#' @export
cal_weightMatrix<-function(trList){
    concept.weight<-c()
    concept.weightDel<-c()
    concept.name<-c()
    print("Calculating weight matrix")
    for(i in seq(1,length(file.concept))){
        tmp.line<-unlist(strsplit(file.concept[[i]],"\t"))
        if(length(tmp.line)-2>=5){
            tmp.conceptVec<-tmp.line[seq(3,length(tmp.line))]
            tmp.data<-as.numeric(unlist(strsplit(cal_conceptWeight(tmp.conceptVec,trList),"_")))
            concept.weight<-append(concept.weight,tmp.data[1]/tmp.data[2])
            concept.weightDel<-append(concept.weightDel,(tmp.data[1]-1)/tmp.data[2])
            concept.name<-append(concept.name,tmp.line[1])
        }
        if(i %% 10000==0){
            print(paste("Processed ",i," lines",sep=""))
        }
    }
    concept.tmpInfo<-data.frame(concept.name,concept.weight,concept.weightDel)
    return(concept.tmpInfo)
}
    
#' @title calculate uniConSig using pre-calculated files
#' @param trList A vector containing all the training gene ids.
#' @param preCal Default using a small subset of the whole pre-calculation file. To do full calculation, see the examples
#' @description calculate the functional associations between a single gene and a training gene set stored in trList.
#' @return A matrix with 3 columns: gene ID, gene Name, and uniConSig score
#' @examples
#' #For customized training list:
#' trList.call<-1:10 #read customized training list
#' result<-cal_uniConSig(trList.call) #calculate uniConSig
#' #To do full calculation:
#' #preCal.local<-get_data_uniConSigPreCal()
#' #load(preCal.local)
#' #result<-cal_uniConSig(trList.call,preCal=preCal.data.all)
#' #head(result) #examine the result
#' @export 
cal_uniConSig<-function(trList,preCal=preCal.data){
    gene.score<-c()
    gene.id<-c()
    concept.info<-cal_weightMatrix(trList)
    print("Calculating uniConSig scores")
    for(j in seq(1,length(preCal))){
        tmp.line<-unlist(strsplit(preCal[[j]],"\t"))
        concept.name<-c()
        concept.epsilon<-c()
        for(i in seq(3,length(tmp.line))){
            tmp.info<-unlist(strsplit(tmp.line[i],"_"))
            if(length(tmp.info)>2){
                concept.name<-append(concept.name,paste(tmp.info[seq(1,(length(tmp.info)-1))],collapse="_"))
            }else{
                concept.name<-append(concept.name,tmp.info[1])
            }
            concept.epsilon<-append(concept.epsilon,tmp.info[length(tmp.info)])
        }
        concept.toCal<-data.frame(concept.name,concept.epsilon)
        tmp.epsilon<- concept.info[concept.info$concept.name %in% concept.toCal$concept.name,]
        tmp.matrix<-as.matrix(merge(concept.toCal,tmp.epsilon,by="concept.name"))
        if(tmp.line[1] %in% trList){
            tmp.vec<-as.numeric(tmp.matrix[,4])/as.numeric(tmp.matrix[,2])
        }else{
            tmp.vec<-as.numeric(tmp.matrix[,3])/as.numeric(tmp.matrix[,2])
        }
        tmp.sum<-sum(tmp.vec)
        tmp.score<-tmp.sum/sqrt(as.numeric(tmp.line[2]))
        gene.score<-append(gene.score,tmp.score)
        gene.id<-append(gene.id,tmp.line[1])
        #print(paste(tmp.line[1],tmp.score,sep="\t"))
        j=j+1
        if(j %% 10000==0){
            print(paste("Processed ",j," genes",sep=""))
        }
    }
    result<-data.frame(gene.id,gene.score)
    result<-merge(gene.anno,result,by="gene.id")
    print("Done.Sorting,and scaling...")
    result<-result[order(result$gene.score,decreasing=TRUE),]
    result[result$gene.score=="NaN",3]<-0
    #result<-result[!result$gene.score=="NaN",]
    result$gene.score<-result$gene.score/max(result$gene.score)
    print("Done.")
    return(result)
}

#' @title Fetch pre-calculated data from github for uniConSig calculation
#' @param file.data A string which specifies the name of the data file to be fetched from github. Default is preCal.data.all.RData
#' @param dir.local A string which specifies the name of the directory where the data file will be saved. Need "/" at the end of the string.
#' @description Download the pre-calculated data from github. Returns the name of the downloaded file with the path to it.
#' @return The name of the downloaded file with the path to it.
#' @examples
#' #For the data required in uniConSig calculation:
#' #preCal.local<-get_data_uniConSigPreCal()
#' #load(preCal.local)
#' @export
get_data_uniConSigPreCal<-function(file.data="preCal.data.all.RData",dir.local="./data.uniConSigPreCal/"){
    dir.remote<-"https://github.com/wangxlab/data_uniConSigPreCal/raw/master/"
    if(!dir.exists(dir.local)){
      dir.create(dir.local)
      print(paste("Directory '",dir.local,"' doesn't exist, created new one.",sep=""))
    }
    url.call<-paste(dir.remote,file.data,sep="")
    file.local<-paste(dir.local,file.data,sep="")
    download.file(url.call,file.local)
    return(file.local)
}

##############################################################CSEA#################################################################

#' @title Calculate cumulative normalized weights
#' @param weight A vector of weights. The name of each value in the vector is it's Entrez gene ID
#' @param posList A list of list containing pathways of genes. Each sublist is a pathway (or a concept), containing Entrez gene IDs
#' @description calculate the cumulative normalized weights for a sorted gene list with weights and an on-set gene set
#' @return A vector of cumulative normalized weights
#' @examples 
#' weight.call<-c(0.1,0.2,0.5,0.7,0.9,1)
#' names(weight.call)<-1:6
#' pathway.my<-list()
#' pathway.my[[1]]<-1:10
#' pathway.my[[2]]<-5:15
#' cumuNorWeight<-cal_cumNorWeight(weight.call,pathway.my)
#' @export
cal_cumNorWeight<-function(weight,posList){
    weight<-sort(weight,decreasing=TRUE)
    onIndice<-which(names(weight) %in% posList)
    onList<-weight[onIndice]
    sumUi<-sum(abs(onList))
    stepDown<-1/(length(weight)-length(onList))
    norWeight<-rep(-stepDown,length(weight))
    names(norWeight)<-names(weight)
    norWeight[onIndice]<-abs(weight[onIndice]/sumUi)
    cumNorWeight<-cumsum(norWeight)
    return(cumNorWeight)
}

#' @title Find the highest absolute value of weights
#' @param cumWeight A vector of cumulative normalized weights. The name of each value in the vector is it's Entrez gene ID
#' @description Find the highest running sum of the weighted random walk k-s test to be the enrichment score (ES)
#' @return A vector of cumulative normalized weights
#' @examples 
#' #For a "cumuNorWeight" calculated from cal_cumNorWeight:
#' weight.call<-c(0.1,0.2,0.5,0.7,0.9,1)
#' names(weight.call)<-1:6
#' pathway.my<-list()
#' pathway.my[[1]]<-1:10
#' pathway.my[[2]]<-5:15
#' cumuNorWeight<-cal_cumNorWeight(weight.call,pathway.my)
#' ES<-findES(cumuNorWeight)
#' @export
findES<-function(cumWeight){
    maxES<-max(abs(cumWeight))
    if(length(cumWeight[which(cumWeight==maxES)])==0){
        ES<-cumWeight[which(cumWeight==-maxES)]
    }else{
        ES<-cumWeight[which(cumWeight==maxES)]
    }
    return(ES)
}

#' @title Calculate the permutations of random on-set genes
#' @param weights A vector of weights. The name of each value in the vector is it's Entrez gene ID
#' @param numOnList An positive integer, indicating the number of genes to be randomly selected
#' @param nPermu Number of permutations, default 1000
#' @description For a certain weighted gene list(weights), randomly select a certain number(numOnList) of genes to be on-set genes, and calculate ES
#' @return A numeric vector of enrichment scores of randomly selected genes
#' @examples 
#' weight.call<-c(0.1,0.2,0.5,0.7,0.9,1)
#' names(weight.call)<-1:6
#' test.tmp<-perm_weightedKS_ofCSEA(weight.call,6)
#' @export
perm_weightedKS_ofCSEA<-function(weights,numOnList,nPermu=1000){
    permu<-c()
    while(length(permu)<nPermu){
        tmpFusList<-sample(names(weights),numOnList)
        permCumNorWeight<-cal_cumNorWeight(weights,tmpFusList)
        permES<-findES(permCumNorWeight)
        permu<-append(permu,permES)
    }
    permu<-as.numeric(permu)
    return(permu)
}

#' @title contruct permutation list
#' @param weights.cp A vector of weights. The name of each value in the vector is it's Entrez gene ID
#' @param posList A list of list, containing pathway genes. See example pathway.hallmark and pathway.c2cp
#' @param nPermu Number of permutations. Default 1000
#' @description contruct permutation list to avoid repeated calculations of same number of randomly selected genes
#' @return A list of list, containing all the enrichment scores of each number of randomly selected genes.
#' @examples 
#' weight.call<-c(0.1,0.2,0.5,0.7,0.9,1)
#' names(weight.call)<-1:6
#' pathway.my<-list()
#' pathway.my[[1]]<-1:10
#' pathway.my[[2]]<-5:15
#' myPermu.call<-construct_permuList(weight.call,pathway.my)
#' @export
construct_permuList<-function(weights.cp,posList,nPermu=1000){
    myPermu<-list()
    n<-0
    for(k in seq(1,length(posList))){
        numOnList.call<-length(weights.cp[which(names(weights.cp) %in% posList[[k]])])
        if(numOnList.call<6){
        }else{
            tmp.index<-which(names(myPermu)==numOnList.call)
            if(length(tmp.index)==0){
                n<-n+1
                #print(paste("No preCal permutation data, calculate on-set gene number ",numOnList.call,sep=""))
                myPermu[[n]]<-perm_weightedKS_ofCSEA(weights.cp,numOnList.call,nPermu = 1000)
                names(myPermu)[n]<-numOnList.call
                #print("Done")
            }
        }
        if(k %% 100==0){
            print(paste("Finished ",k," posLists",sep=""))
        }
    }
    return(myPermu)
}

#' @title Calculate weighted random walk k-s test for CSEA
#' @param weights A vector of weights. The name of each value in the vector is it's Entrez gene ID
#' @param positiveList A list of list containing pathways of genes. Each sublist is a pathway (or a concept), containing Entrez gene IDs
#' @param myPermu A list of list containing permutations of randomly selected genes for each pathway.
#' @param nPermu The number of permutations used throughout the calculations
#' @description calculate normalized enrichment score based on weighted gene list(weights), on-set gene set(positiveList), and permutation result(myPermu).
#' @return A table of 3 columns, "pathway names","NES","pValue"
#' @examples 
#' ##For a "result" calculated from cal_uniConSig:
#' weight.call<-c(0.1,0.2,0.5,0.7,0.9,1)
#' names(weight.call)<-1:6
#' pathway.my<-list()
#' pathway.my[[1]]<-1:10
#' pathway.my[[2]]<-5:15
#' myPermu.call<-construct_permuList(weight.call,pathway.my)
#' myNES<-weightedKS_ofCSEA(weight.call,pathway.my,myPermu.call)
#' @export
weightedKS_ofCSEA<-function(weights,positiveList,myPermu,nPermu=1000){
    NESResult<-c()
    pVResult<-c()
    for(k in seq(1,length(positiveList))){
        myCumNorWeight<-cal_cumNorWeight(weights,positiveList[[k]])
        ES<-findES(myCumNorWeight)
        numOnList<-length(weights[which(names(weights) %in% positiveList[[k]])])
        if(numOnList<6){
            NES<-"NA"
            names(NES)<-names(positiveList[k])
            NESResult<-append(NESResult,NES)
            pValue<-"NA"
            names(pValue)<-names(positiveList[k])
            pVResult<-append(pVResult,pValue)
        }else{
            tmp.index<-which(names(myPermu)==numOnList)
            if(length(tmp.index)==0){
                print(paste("No preCal permutation data, calculate number ",numOnList,sep=""))
                permu<-perm_weightedKS_ofCSEA(weights = weights,numOnList = numOnList,nPermu = nPermu)
                tmp.num<-length(names(myPermu))+1
                myPermu[[tmp.num]]<-permu
                names(myPermu)[tmp.num]<-numOnList
                print(paste("Adding new permutation data of on-set number:",numOnList,sep=""))
                print("Done")
            }else{
                permu<-myPermu[[tmp.index]]
            }
            NES<-ES/mean(abs(permu))
            names(NES)<-names(positiveList[k])
            NESResult<-append(NESResult,NES)
            pValue<-length(permu[which(permu>ES)])/nPermu
            names(pValue)<-names(positiveList[k])
            if(NES>0){
                pVResult<-append(pVResult,pValue)
            }else{
                pVResult<-append(pVResult,1-pValue)
            }
        }
        if(k %% 10000==0){
            print(paste("Finished ",k," posLists",sep=""))
        }
    }
    myResult<-merge(NESResult,pVResult,by="row.names")
    colnames(myResult)=c("names","NES","pValue")
    return(myResult)
}

#' @title Calculate CSEA
#' @param result.uniConSig The result generated by "cal_uniConSig"
#' @param posList A list of list containing pathways of genes. Each sublist is a pathway (or a concept), containing Entrez gene IDs
#' @param nPermu The number of permutations. Default 1000. Note that this number should be consistent with pre-calculated permutation list.
#' @description transfer the uniConSig result to function "weightedKS_ofCSEA()"
#' @return A table of 3 columns, "pathway names","NES","pValue"
#' @examples 
#' ##For a training list, starting from uniConSig:
#' trList.call<-1:100
#' result<-cal_uniConSig(trList.call)
#' pathway.my<-list()
#' pathway.my[[1]]<-1:10
#' pathway.my[[2]]<-5:15
#' result.csea<-CSEA(result,pathway.my)
#' @export
CSEA<-function(result.uniConSig,posList,nPermu=1000){
    result.uniConSig<-result.uniConSig[!result.uniConSig$gene.score==0,]
    tmp.weight<-result.uniConSig$gene.score
    names(tmp.weight)<-result.uniConSig$gene.id
    print("Calculating permutations...")
    preCal.permutations<-construct_permuList(tmp.weight,posList)
    print("Done.")
    print("Calculating NES")
    myFinal<-weightedKS_ofCSEA(tmp.weight,posList,preCal.permutations,nPermu)
    print("Done")
    myFinal<-myFinal[order(myFinal$NES,decreasing = TRUE),]
    return(myFinal)
}

#Documentation of datasets
#For pathway.c2cp
#' @title Pathway gene sets of C2CP from MSigDB
#' @name pathway.c2cp
#' @description Compiled gene sets from MSigdDB c2cp
#' @docType data
#' @author Xu Chi (original data from MSigDB c2cp) \email{xuc13@@pitt.edu}
#' @keywords datasets
#' @usage data(pathway.c2cp)
#' @format A list of list, elements are entrez gene ids. Names of each sublist is the pathway's name
NULL

#For pathway.hallmark
#' @title Pathway gene sets of HALLMARK from MSigDB
#' @name pathway.hallmark
#' @description Compiled gene sets from MSigDB hallmark
#' @docType data
#' @author Xu Chi (original data from MSigDB c2cp) \email{xuc13@@pitt.edu}
#' @keywords datasets
#' @usage data(pathway.hallmark)
#' @format A list of list, elements are entrez gene ids. Names of each sublist is the pathway's name
NULL

#For trList.cgc
#' @title Cancer gene sets from Cancer Gene Census
#' @name trList.cgc
#' @description A list of cancer genes generated by Cancer Gene Census
#' @docType data
#' @author Xu Chi (original data from Cancer Gene Census) \email{xuc13@@pitt.edu}
#' @keywords datasets
#' @usage data(trList.cgc)
#' @format A list entrez gene ids.
NULL
