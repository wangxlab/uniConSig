### Input File Format ####
###  Input file should have gene IDs or symbols in row, and sample names in column. 

#####################################################################
######################### Package to install ###########################
#install.packages("plyr")
library(plyr)

###########################################################################
##################### User modification ##############################
setwd("Y:/CaGenome/Pipelines/uniConSig/ConvertGeneExpData_to_GMT")
inputFile="ExpressionDataExample.txt"     
outputFile="ExpressionDataGMT_Example.gmt"
scaleFactor=1.4826    # median+(scaleFactor*MAD)
MAD_constant=1      # This is constant for MAD calculation
options(digits=5)   # To apply decimal numbers. 
##################################################################

expData=read.table(inputFile, stringsAsFactors=F, row.names=1, header=T, sep="\t")  #Read input data as data.frame
#expDataGeneSym<-cbind(expData, rownames(expData))

smpByMAD<-data.frame(stringsAsFactors=F)   # Initialize a data.frame
selectSmp<-data.frame(stringsAsFactors=F)  # Initialize a data.frame

#####################################################################
### function to find samples that have gene expression level > median + (scale factor * MAD)  or gene expression level < median + (scale factor * MAD)
madCalculate<-function(x, constantGiven=MAD_constant) {
	#constantGiven=1
	#x<-expData[1,]
	#medianValue<-median(as.numeric(x))
	#madValue<-mad(as.numeric(x), constant=constantGiven, na.rm=T) 
	

	medianValue<-median(as.numeric(x))    # I shoud make x as.numeric
	print(paste("median: ", medianValue, sep=""))
	
	madValue<-mad(as.numeric(x), constant=constantGiven, na.rm=T)          # I shoud make x as.numeric
	print(paste("MAD: ", madValue, sep=""))
	
	highCriteria<-medianValue+(scaleFactor*madValue)
	lowCriteria<-medianValue-(scaleFactor*madValue)
	upSmp<-which (x>highCriteria)
	downSmp<-which(x<lowCriteria)	

	upSmpName<-names(expData)[upSmp]
	downSmpName<-names(expData)[downSmp]
	length(upSmpName)<-max(length(upSmpName), length(downSmpName))
	length(downSmpName)<-max(length(upSmpName), length(downSmpName))	

	selectSmp<-rbind(upSmpName, downSmpName)
	smpByMAD<-rbind(smpByMAD, selectSmp)

	return(smpByMAD)
}

## Calculate MAD in each row of data, and find samples that have gene expression higher or lower than the criteria.
sampleByMAD<-apply(expData, 1, madCalculate)   

gmtFormDraft<-data.frame(stringsAsFactors=F)  # Initialize a data.frame
oneGeneGMT<-data.frame(stringsAsFactors=F)    # Initialize a data.frame

##########################################################################
### From the output of "madCalculate" function, make gmt format
##########################################################################

len_smpByMAD<-length(sampleByMAD)
for (i in seq(len_smpByMAD)) {
	#i=1
	geneSym<-names(sampleByMAD[i])
	geneSym_Up<-paste(geneSym, "_up", sep="")
	geneSym_Down<-paste(geneSym, "_down", sep="")

	sampleUp<-NULL
	sampleDown<-NULL
	if (nrow(sampleByMAD[[i]]==2))	{
		smpUpFound<-!is.na(sampleByMAD[[i]][1,])
		smpDownFound<-!is.na(sampleByMAD[[i]][2,])
		sampleUp<-cbind(geneSym_Up, length(smpUpFound[smpUpFound==TRUE]), as.matrix(sampleByMAD[[i]][1,]))
		sampleDown<-cbind(geneSym_Down, length(smpDownFound[smpDownFound==TRUE]), as.matrix(sampleByMAD[[i]][2,]))
	} else if (nrow(sampleByMAD[[i]])==1) {
		if (rownames(sampleByMAD[[i]])=="upSmpName") {
			smpUpFound<-!is.na(sampleByMAD[[i]][1,])
			sampleUp<-cbind(geneSym_Up, length(smpUpFound[smpUpFound==TRUE]), as.matrix(sampleByMAD[[i]][1,]))
		} else if (rownames(sampleByMAD[[i]])=="downSmpName") {
			smpDownFound<-!is.na(sampleByMAD[[i]][1,])
			sampleDown<-cbind(geneSym_Down, length(smpDownFound[smpDownFound==TRUE]), as.matrix(sampleByMAD[[i]][1,]))
		}
	} 

	colnames(sampleUp)<-NULL
	colnames(sampleDown)<-NULL

	print(paste("sampleUp: ", sampleUp)) 
	print(paste("sampleDown: ", sampleDown))
	
	oneGeneGMT<-rbind(sampleUp, sampleDown)
	oneGeneGMT_df<-as.data.frame(oneGeneGMT)	

	gmtFormDraft<-rbind.fill(gmtFormDraft, oneGeneGMT_df)	
}

gmtForm<-gmtFormDraft[-which(is.na(gmtFormDraft[,3])),]   # how do I delete row in data.frame: myData <- myData[-c(2, 4, 6), ]

#######################################################################################
## If the sample names starts with numeric, so it generated "X" in the beginning of the samples names, 
## I may want to remove "X" in the sample names
#######################################################################################
removeX<-function (dataMatrix) {
	#print(dataMatrix)
	smpNameRemX<-gsub("^X", "", dataMatrix)
	#print(smpNameRemX)
	return(smpNameRemX)
}
gmtFormRemX<-apply(gmtForm[,3:ncol(gmtForm)], 2, removeX)
gmtFormGeneSymSize<-cbind(gmtForm[,1:2], gmtFormRemX)
###########################################

gmtFormGeneSymSizeChar <- sapply(gmtFormGeneSymSize, as.character)   # If values are 'factor', I can't change NA to blank
gmtFormGeneSymSizeChar[is.na(gmtFormGeneSymSizeChar)] <- ""   # change 'na' to empty 

write.table (gmtFormGeneSymSizeChar, outputFile, quote=F, sep="\t", col.names=F, row.names=F)

###############################################################################
### Personal practice 
###############################################################################


#gmtFormChar <- sapply(gmtForm, as.character) # if values are `factor`,
#gmtFormChar[is.na(gmtFormChar)] <- ""   # change 'na' to empty 
#gmtFormChar[is.na(gmtFormChar)] <- ""   # change 'na' to empty in a matrix. 
#write.table (gmtFormFinal, outputFile, quote=F, sep="\t", col.names=F, row.names=F)
#gmtFormDF<-as.data.frame(gmtFormChar)   # make gmtFormChar data.frame
#write.table (gmtFormFinal, outputFile, quote=F, sep="\t", col.names=F, row.names=F)

