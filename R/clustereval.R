### function cluster_performance() computes external cluster performance measures: Purity, V-Measure, Normalised Mutual Information and Adjusted Rand Index
## arguments:
# assigned = vector of cluster assignments
# labels = true class labels (same length as assigned)
# beta = weight parameter used in calculation of V-measure. Higher values apply higher weight to homogeneity over completeness.

cluster_performance = function(assigned, labels, beta = 1){
	n <- length(labels)
	T <- table(assigned, labels)
	RS <- rowSums(T)
	CS <- colSums(T)

	## V-measure

	CK <- - sum(apply(T, 1, function(x) return(sum(x[which(x>0)]*log(x[which(x>0)]/sum(x))))))/n
	KC <- - sum(apply(T, 2, function(x) return(sum(x[which(x>0)]*log(x[which(x>0)]/sum(x))))))/n
	K <- - sum(apply(T, 1, function(x) return(sum(x)*log(sum(x)/n))))/n
	C <- - sum(apply(T, 2, function(x) return(sum(x)*log(sum(x)/n))))/n
	if(C!=0){
		h <- 1 - CK/C
	}
	else{
		h <- 0
	}
	if(K!=0){
		c <- 1 - KC/K
	}
	else{
		c <- 0
	}
	if(h==0 && c==0) v.measure <- 0
	else v.measure <- (1+beta)*h*c/(beta*h+c)

	## Purity

	purity <- sum(apply(T, 1, function(x) return(max(x))))/n

	## Adjusted Rand Index

	O <- sum(sapply(T, function(t) choose(t, 2)))
	E <- (sum(sapply(RS, function(t) choose(t, 2)))*sum(sapply(CS, function(t) choose(t, 2))))/choose(n, 2)
	M <- (sum(sapply(RS, function(t) choose(t, 2))) + sum(sapply(CS, function(t) choose(t, 2))))/2
	adj.rand <- (O-E)/(M-E)


	## Normalised Mutual Information

	prod <- RS%*%t(CS)
	Tp <- T
	Tp[which(T==0)] <- 1e-10
	IXY <- sum(T*log(Tp*n/prod))
	HX <- sum(RS*log(RS/n))
	HY <- sum(CS*log(CS/n))
	NMI <- IXY/sqrt(HX*HY)

	c(adj.rand = adj.rand, purity = purity, v.measure = v.measure, nmi = NMI)
}



### function success_ratio() computes the success ratio (Pavlidis et al. 2016) of a binary partition
## arguments:
# assigned = vector of cluster assignments (takes at most two distinct values)
# labels = true class labels (same length as assigned)



success_ratio = function(assigned, labels){
  n = length(labels)

  #### Find the sizes of the overlaps between assigned and labels
  T = table(assigned, labels)

  #### determine which of the 2 aggregated classes to assign the classes
  split = apply(T, 2, which.max)

  #### if an aggregated class is empty, exit. failure to split any classes
  if(length(unique(split))==1){
    return(0)
  }

  #### calculate success and error
  success = min(sum(T[1,split==1]), sum(T[2,split==2]))
  error = sum(apply(T, 2, min))

  success/(success+error)
}
