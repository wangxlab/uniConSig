library("dplyr")
# The Broad Institute SOFTWARE COPYRIGHT NOTICE AGREEMENT This software and its
# documentation are copyright 2003 by the Broad Institute/Massachusetts Institute
# of Technology.  All rights are reserved.  This software is supplied without any
# warranty or guaranteed support whatsoever. Neither the Broad Institute nor MIT
# can be responsible for its use, misuse, or functionality.

# G S E A -- Gene Set Enrichment Analysis
# Auxiliary functions and definitions

GSEA.ReadClsFile <- function(file = "NULL") {
 # Reads a class vector CLS file and defines phenotype and class labels vectors
 # for the samples in a gene expression file (RES or GCT format)
 
 cls.cont <- readLines(file)
 num.lines <- length(cls.cont)
 cls.cont[[3]] <- gsub("\\t", " ", cls.cont[[3]])  #Converts any tabs to spaces
 class.list <- unlist(strsplit(cls.cont[[3]], " "))  #Splits CLS on spaces
 s <- length(class.list)
 t <- table(class.list)
 l <- length(t)
 phen <- vector(length = l, mode = "character")
 phen.label <- vector(length = l, mode = "numeric")
 class.v <- vector(length = s, mode = "numeric")
 for (i in 1:l) {
  phen[i] <- noquote(names(t)[i])
  phen.label[i] <- i - 1
 }
 for (i in 1:s) {
  for (j in 1:l) {
   if (class.list[i] == phen[j]) {
    class.v[i] <- phen.label[j]
   }
  }
 }
 return(list(phen = phen, class.v = class.v))
}



GSEA.GeneRanking <- function(A, class.labels, gene.labels, nperm, permutation.type = 0, 
 sigma.correction = "GeneCluster", fraction = 1, replace = F, reverse.sign=F, rank.metric) {
 
 A <- A + 1e-08
 
 N <- length(A[, 1])
 Ns <- length(A[1, ])
 
 subset.mask <- matrix(0, nrow = Ns, ncol = nperm)
 reshuffled.class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
 reshuffled.class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
 class.labels1 <- matrix(0, nrow = Ns, ncol = nperm)
 class.labels2 <- matrix(0, nrow = Ns, ncol = nperm)
 
 order.matrix <- matrix(0, nrow = N, ncol = nperm)
 obs.order.matrix <- matrix(0, nrow = N, ncol = nperm)
 rnk.matrix <- matrix(0, nrow = N, ncol = nperm)
 obs.rnk.matrix <- matrix(0, nrow = N, ncol = nperm)
 
 obs.gene.labels <- vector(length = N, mode = "character")
 obs.gene.descs <- vector(length = N, mode = "character")
 obs.gene.symbols <- vector(length = N, mode = "character")
 
 M1 <- matrix(0, nrow = N, ncol = nperm)
 M2 <- matrix(0, nrow = N, ncol = nperm)
 S1 <- matrix(0, nrow = N, ncol = nperm)
 S2 <- matrix(0, nrow = N, ncol = nperm)
 
 gc()
 
 C <- split(class.labels, class.labels)
 class1.size <- length(C[[1]])
 class2.size <- length(C[[2]])
 class1.index <- seq(1, class1.size, 1)
 class2.index <- seq(class1.size + 1, class1.size + class2.size, 1)
 
 for (r in 1:nperm) {
  class1.subset <- sample(class1.index, size = ceiling(class1.size * fraction), 
   replace = replace)
  class2.subset <- sample(class2.index, size = ceiling(class2.size * fraction), 
   replace = replace)
  class1.subset.size <- length(class1.subset)
  class2.subset.size <- length(class2.subset)
  subset.class1 <- rep(0, class1.size)
  for (i in 1:class1.size) {
   if (is.element(class1.index[i], class1.subset)) {
    subset.class1[i] <- 1
   }
  }
  subset.class2 <- rep(0, class2.size)
  for (i in 1:class2.size) {
   if (is.element(class2.index[i], class2.subset)) {
    subset.class2[i] <- 1
   }
  }
  subset.mask[, r] <- as.numeric(c(subset.class1, subset.class2))
  fraction.class1 <- class1.size/Ns
  fraction.class2 <- class2.size/Ns
  
  if (permutation.type == 0) {
   # random (unbalanced) permutation
   full.subset <- c(class1.subset, class2.subset)
   label1.subset <- sample(full.subset, size = Ns * fraction.class1)
   reshuffled.class.labels1[, r] <- rep(0, Ns)
   reshuffled.class.labels2[, r] <- rep(0, Ns)
   class.labels1[, r] <- rep(0, Ns)
   class.labels2[, r] <- rep(0, Ns)
   for (i in 1:Ns) {
    m1 <- sum(!is.na(match(label1.subset, i)))
    m2 <- sum(!is.na(match(full.subset, i)))
    reshuffled.class.labels1[i, r] <- m1
    reshuffled.class.labels2[i, r] <- m2 - m1
    if (i <= class1.size) {
      class.labels1[i, r] <- m2
      class.labels2[i, r] <- 0
    } else {
      class.labels1[i, r] <- 0
      class.labels2[i, r] <- m2
    }
   }
  } else if (permutation.type == 1) {
   # proportional (balanced) permutation
   
   class1.label1.subset <- sample(class1.subset, size = ceiling(class1.subset.size * 
    fraction.class1))
   class2.label1.subset <- sample(class2.subset, size = floor(class2.subset.size * 
    fraction.class1))
   reshuffled.class.labels1[, r] <- rep(0, Ns)
   reshuffled.class.labels2[, r] <- rep(0, Ns)
   class.labels1[, r] <- rep(0, Ns)
   class.labels2[, r] <- rep(0, Ns)
   for (i in 1:Ns) {
    if (i <= class1.size) {
      m1 <- sum(!is.na(match(class1.label1.subset, i)))
      m2 <- sum(!is.na(match(class1.subset, i)))
      reshuffled.class.labels1[i, r] <- m1
      reshuffled.class.labels2[i, r] <- m2 - m1
      class.labels1[i, r] <- m2
      class.labels2[i, r] <- 0
    } else {
      m1 <- sum(!is.na(match(class2.label1.subset, i)))
      m2 <- sum(!is.na(match(class2.subset, i)))
      reshuffled.class.labels1[i, r] <- m1
      reshuffled.class.labels2[i, r] <- m2 - m1
      class.labels1[i, r] <- 0
      class.labels2[i, r] <- m2
    }
   }
  }
 }
 
 # compute S2N for the random permutation matrix
 if (rank.metric == "S2N") {
  P <- reshuffled.class.labels1 * subset.mask
  n1 <- sum(P[, 1])
  M1 <- A %*% P
  M1 <- M1/n1
  gc()
  A2 <- A * A
  S1 <- A2 %*% P
  S1 <- S1/n1 - M1 * M1
  S1 <- sqrt(abs((n1/(n1 - 1)) * S1))
  gc()
  P <- reshuffled.class.labels2 * subset.mask
  n2 <- sum(P[, 1])
  M2 <- A %*% P
  M2 <- M2/n2
  gc()
  A2 <- A * A
  S2 <- A2 %*% P
  S2 <- S2/n2 - M2 * M2
  S2 <- sqrt(abs((n2/(n2 - 1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {
   # small sigma 'fix' as used in GeneCluster
   S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
   S2 <- ifelse(S2 == 0, 0.2, S2)
   S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
   S1 <- ifelse(S1 == 0, 0.2, S1)
   gc()
  }
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()
  
  rnk.matrix <- M1/S1
  
  if (reverse.sign == T) {
   rnk.matrix <- -rnk.matrix
  }
  gc()
  
  for (r in 1:nperm) {
   order.matrix[, r] <- order(rnk.matrix[, r], decreasing = T)
  }
  
  # compute S2N for the 'observed' permutation matrix
  
  P <- class.labels1 * subset.mask
  n1 <- sum(P[, 1])
  M1 <- A %*% P
  M1 <- M1/n1
  gc()
  A2 <- A * A
  S1 <- A2 %*% P
  S1 <- S1/n1 - M1 * M1
  S1 <- sqrt(abs((n1/(n1 - 1)) * S1))
  gc()
  P <- class.labels2 * subset.mask
  n2 <- sum(P[, 1])
  M2 <- A %*% P
  M2 <- M2/n2
  gc()
  A2 <- A * A
  S2 <- A2 %*% P
  S2 <- S2/n2 - M2 * M2
  S2 <- sqrt(abs((n2/(n2 - 1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {
   # small sigma 'fix' as used in GeneCluster
   S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
   S2 <- ifelse(S2 == 0, 0.2, S2)
   S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
   S1 <- ifelse(S1 == 0, 0.2, S1)
   gc()
  }
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- S1 + S2
  rm(S2)
  gc()
  
  obs.rnk.matrix <- M1/S1
  gc()
 } else if (rank.metric == "ttest") {
 # compute TTest for the random permutation matrix
  P <- reshuffled.class.labels1 * subset.mask
  n1 <- sum(P[, 1])
  M1 <- A %*% P
  M1 <- M1/n1
  gc()
  A2 <- A * A
  S1 <- A2 %*% P
  S1 <- S1/n1 - M1 * M1
  S1 <- sqrt(abs((n1/(n1 - 1)) * S1))
  gc()
  P <- reshuffled.class.labels2 * subset.mask
  n2 <- sum(P[, 1])
  M2 <- A %*% P
  M2 <- M2/n2
  gc()
  A2 <- A * A
  S2 <- A2 %*% P
  S2 <- S2/n2 - M2 * M2
  S2 <- sqrt(abs((n2/(n2 - 1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {
   # small sigma 'fix' as used in GeneCluster
   S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
   S2 <- ifelse(S2 == 0, 0.2, S2)
   S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
   S1 <- ifelse(S1 == 0, 0.2, S1)
   gc()
  }
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- (S1^2)/class1.size
  S2 <- (S2^2)/class2.size
  S1 <- S1 + S2
  S1 <- sqrt(S1)
  rm(S2)
  gc()
  
  rnk.matrix <- M1/S1
  
  if (reverse.sign == T) {
   rnk.matrix <- -rnk.matrix
  }
  gc()
  
  for (r in 1:nperm) {
   order.matrix[, r] <- order(rnk.matrix[, r], decreasing = T)
  }
  
  # compute TTest for the 'observed' permutation matrix
  
  P <- class.labels1 * subset.mask
  n1 <- sum(P[, 1])
  M1 <- A %*% P
  M1 <- M1/n1
  gc()
  A2 <- A * A
  S1 <- A2 %*% P
  S1 <- S1/n1 - M1 * M1
  S1 <- sqrt(abs((n1/(n1 - 1)) * S1))
  gc()
  P <- class.labels2 * subset.mask
  n2 <- sum(P[, 1])
  M2 <- A %*% P
  M2 <- M2/n2
  gc()
  A2 <- A * A
  S2 <- A2 %*% P
  S2 <- S2/n2 - M2 * M2
  S2 <- sqrt(abs((n2/(n2 - 1)) * S2))
  rm(P)
  rm(A2)
  gc()
  
  if (sigma.correction == "GeneCluster") {
   # small sigma 'fix' as used in GeneCluster
   S2 <- ifelse(0.2 * abs(M2) < S2, S2, 0.2 * abs(M2))
   S2 <- ifelse(S2 == 0, 0.2, S2)
   S1 <- ifelse(0.2 * abs(M1) < S1, S1, 0.2 * abs(M1))
   S1 <- ifelse(S1 == 0, 0.2, S1)
   gc()
  }
  
  M1 <- M1 - M2
  rm(M2)
  gc()
  S1 <- (S1^2)/class1.size
  S2 <- (S2^2)/class2.size
  S1 <- S1 + S2
  S1 <- sqrt(S1)
  rm(S2)
  gc()
  
  obs.rnk.matrix <- M1/S1
  gc()
 }
 
 if (reverse.sign == T) {
  obs.rnk.matrix <- -obs.rnk.matrix
 }
 
 for (r in 1:nperm) {
  obs.order.matrix[, r] <- order(obs.rnk.matrix[, r], decreasing = T)
 }
 
 return(list(rnk.matrix = rnk.matrix, obs.rnk.matrix = obs.rnk.matrix, order.matrix = order.matrix, 
  obs.order.matrix = obs.order.matrix))
}
GSEA_SignalToNoise  <- function(input.ds, input.cls, nperm=1000) { 
	reverse.sign = F
	preproc.type = 0
	random.seed=as.integer(as.POSIXct(Sys.time())) 
	perm.type = 0
	fraction = 1
	replace = F
	rank.metric="S2N"

	print(" *** Running Gene Set Enrichment Analysis...")
	# Start of GSEA methodology
	# Read input data matrix
 
	set.seed(seed = random.seed, kind = NULL)
	adjust.param <- 0.5
 
	time1 <- proc.time()

	#dataset <- GSEA.Gct2Frame(filename = input.ds)
	dataset <-  read.table(input.ds,sep="\t",comment.char="",quote="",stringsAsFactors=FALSE,fill=TRUE,header=F)
	dataset <- dataset[-1, ]
	dataset <- dataset[-1, ]
	colnames(dataset) <- dataset[1, ]
	dataset <- dataset[-1, ]  # [1] 38175    42
 
	rownames(dataset) <- dataset[, 1]
	dataset <- dataset[, -1]
	dataset <- dataset[, -1]
	#print(paste("dataset: ", dim(dataset)))  # "dataset:  15645" "dataset:  40"
	#print(dataset[1:5,1:5])
       
 	#print ("datast done")
	gene.labels <- row.names(dataset)
	sample.names <- colnames(dataset[2:length(colnames(dataset))])
	A <- data.matrix(dataset)
	print(paste("A:  ", dim(A)))
  	print(A[1:5,1:5])
	cols <- length(A[1, ])
	rows <- length(A[, 1])
  
	if (is.list(input.cls)) {
		CLS <- input.cls
	} else {
		CLS <- GSEA.ReadClsFile(file = input.cls)
 	}
	class.labels <- CLS$class.v
	class.phen <- CLS$phen
  
	if (reverse.sign == T) {
		phen1 <- class.phen[2]
		phen2 <- class.phen[1]
	} else {
		phen1 <- class.phen[1]
		phen2 <- class.phen[2]
 	}
  
	# sort samples according to phenotype
	col.index <- order(class.labels, decreasing = F)
	class.labels <- class.labels[col.index]
	sample.names <- sample.names[col.index]
	for (j in 1:rows) {
		A[j, ] <- A[j, col.index]
	}
	names(A) <- sample.names
 
	N <- length(A[, 1])
	Ns <- length(A[1, ])
 
	print(c("Number of genes:", N))
	print(c("Number of samples:", Ns))
 
	time2 <- proc.time()
 
	# GSEA methodology 
	# Compute observed and random permutation gene rankings
 
 	obs.rnk <- vector(length = N, mode = "numeric")
	correl.matrix <- matrix(nrow = N, ncol = nperm)
	obs.correl.matrix <- matrix(nrow = N, ncol = nperm)
	order.matrix <- matrix(nrow = N, ncol = nperm)
	obs.order.matrix <- matrix(nrow = N, ncol = nperm)
 
	nperm.per.call <- 100
	n.groups <- nperm%/%nperm.per.call
	n.rem <- nperm%%nperm.per.call
	n.perms <- c(rep(nperm.per.call, n.groups), n.rem)
	n.ends <- cumsum(n.perms)
	n.starts <- n.ends - n.perms + 1
 
	if (n.rem == 0) {
		n.tot <- n.groups
	} else {
		n.tot <- n.groups + 1
	}

	for (nk in 1:n.tot) {
		call.nperm <- n.perms[nk]
   		print(paste("Computing ranked list for actual and permuted phenotypes.......permutations: ", 
		n.starts[nk], "--", n.ends[nk], sep = " "))
   
		O <- GSEA.GeneRanking(A, class.labels, gene.labels, call.nperm, permutation.type = perm.type, 
			sigma.correction = "GeneCluster", fraction = fraction, replace = replace, 
			reverse.sign = reverse.sign, rank.metric)
		gc()
   
  		order.matrix[, n.starts[nk]:n.ends[nk]] <- O$order.matrix
		obs.order.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.order.matrix
		correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$rnk.matrix
		obs.correl.matrix[, n.starts[nk]:n.ends[nk]] <- O$obs.rnk.matrix	
		geneName<-rownames(O$obs.rnk.matrix)	
		rm(O)
	}
	obs.rnk <- apply(obs.correl.matrix, 1, median)  # using median to assign enrichment scores
	names(obs.rnk)<-geneName
	obs.rnk <- obs.rnk*-1
	obs.rnk <- sort(obs.rnk, decreasing = T)
	return(obs.rnk)
 
} # end of GSEA_SignalToNoise function