library(data.table)
library(dplyr)
library(GSEA)

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
	

func_makegmt<-function(topPathway.del) {
   #listData<-topPathway.del 
  for (listIndex in seq(1:length(topPathway.del))) {
      geneSigName<-names(topPathway.del[listIndex])
      website<-paste("http://www.broadinstitute.org/gsea/msigdb/cards/", geneSigName,sep="")
      geneSymb<-unname(unlist(topPathway.del[listIndex]))
      gmtRow<-c(geneSigName,website,geneSymb)
      gmtRowMtrx<-t(as.matrix(gmtRow))
      write.table(gmtRowMtrx, "topPathway.del.gmt",append=T,col.names=F,row.names=F,quote=F,sep="\t")  # strange. row.names=T doesn't work.
  }
}

func_findGSEAout2<-function(GSEAout_onefile) {
	## Find GSEA output files	
	#GSEA.resFile<-list.files(path=output.directory, pattern=glob2rx("*.SUMMARY.RESULTS.REPORT.*.txt"))
	#print(paste(output.directory, GSEA.resFile[grep(subtype,GSEA.resFile)], sep=""))
	#GSEA.res<-fread(paste(output.directory,"/", GSEA.resFile[grep(subtype,GSEA.resFile)], sep=""),stringsAsFactors=F,header=T,check.names=F)
	GSEA.result<-GSEAout_onefile %>% dplyr::select(GS,NES,paste("NOM","p-val",sep=" "),paste("FDR","q-val",sep=" ")) %>% rename(Compare.List=GS,NES=NES,pValue=paste("NOM","p-val",sep=" "),qValue=paste("FDR","q-val",sep=" "))
	#GSEA.result$Compare.List<-as.character(GSEA.result$Compare.List)
	GSEA.result$Compare.List<-as.character(levels(GSEA.result$Compare.List))[GSEA.result$Compare.List]
	GSEA.result$NES<-as.numeric(levels(GSEA.result$NES))[GSEA.result$NES]
	GSEA.result$pValue<-as.numeric(levels(GSEA.result$pValue))[GSEA.result$pValue]
	GSEA.result$qValue<-as.numeric(levels(GSEA.result$qValue))[GSEA.result$qValue]
	return(GSEA.result)
}

#disambiguate pathway result
disambiguation_RGSEA<-function(GSEA.result,target.score,gs.db,upPathway,topn=30,p.cut=0.01){ #upPathways TRUE or FALSE
	compare.list<-read_gmt(gs.db,min=10)
	topPathway=GSEA.result$Compare.List[1:topn]  # length(topPathway): 30   ## GSEA.result should not be factor. If it is factor, topPathway is not correctly extracted
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
	target.score.DF<-data.frame(matrix(c(names(target.score),target.score), ncol=2), stringsAsFactors=FALSE ) 
	target.score.DF[,2]<-as.numeric(target.score.DF[,2])
	write.table(target.score.DF, "temp_target.score.txt",col.names=F,row.names=F,sep="\t",quote=F)
	
	## Convert topPathway.del (list) to .gmt and save the file. 	
	func_makegmt(topPathway.del)   # store .gmt file. The file name is "topPathway.del.gmt" # 1,740 pathways
  
	## Run R-GSEA_preranked 
	prerank.output.directory<-"./R-GSEAprerankedOut"
	dir.create(prerank.output.directory)

	## Run GSEA: it takes about 10 minutes
	GSEAoutprerank <- GSEA(   ## Input/Output Files :-------------------------------------------
		input.ds="temp_target.score.txt",	# Input gene expression dataset file in GCT format
		gs.db = "topPathway.del.gmt",		# Gene set database in GMT format
		output.directory=prerank.output.directory,			# It needs "/" at the end. Directory where to store output and results (default: "")
		doc.string="dismb.prerank",   									# Documentation string used as a prefix to name result files (default: "GSEA.analysis")
		reshuffling.type="gene.labels",
		collapse.mode="median",
		gsea.type="preranked"
	)
	unlink("topPathway.del.gmt")
	unlink("temp_target.score.txt")

	if(upPathway) {
		GSEA.del<-func_findGSEAout2(GSEAoutprerank[[1]])   # 801 4    1 4
	} else {
		GSEA.del<-func_findGSEAout2(GSEAoutprerank[[2]])    # 17 4    775 4
	}
	
	tmp.bin <- matrix(0, nrow = length(topPathway), ncol = length(topPathway))
	dimnames(tmp.bin) <- list(topPathway, topPathway)
	ij=which(tmp.bin==0, arr.ind = T)
	for (k in 1:nrow(ij)){
		delname=paste(rownames(tmp.bin)[ij[k,"row"]],colnames(tmp.bin)[ij[k,"col"]],sep=":")
	    if (delname %in% GSEA.del$Compare.List){      # I can't find any delname in GSEA.del$Compare.List ?
			tmp.bin[ij[k,"row"],ij[k,"col"]]=as.numeric(GSEA.del$pValue[GSEA.del$Compare.List==delname])
		}else{
			tmp.bin[ij[k,"row"],ij[k,"col"]]=1
		}
	}

	index=which(tmp.bin>0.01,arr.ind = TRUE)   
	remove=c()
	for (x in 1:nrow(index)){
		rowname=rownames(tmp.bin)[index[x,"row"]]
		colname=colnames(tmp.bin)[index[x,"col"]]
		if (rowname==colname | rowname %in%remove |colname %in% remove ){
		next
		}
		remove=c(remove,rowname)
	}
	filter.result=GSEA.result[1:topn,]
	filter.result=filter.result[which(!filter.result$Compare.List %in% remove),]

	unlink(prerank.output.directory, recursive=T)
	return(list(pathway.filtered=filter.result,pathway.todel=remove,disambiguate.matrix=tmp.bin))
}
