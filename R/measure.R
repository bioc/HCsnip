measure <-
function(parti, dis, X = NULL, method = "g2", maxmiss = 30, ...) 
{
   if (length(parti) != nrow(as.matrix(dis))) 
      stop("input arguments dimension mismatch")
   if (!isSymmetric(as.matrix(dis))) 
      stop("'dis' should be a square distance matrix or 'dist' class object")
   if (method == "c.index") {  #small is good 
      c.index <- c()
      dis <- as.matrix(dis)
      for(i in unique(parti)) {
         s <- subset(dis, (parti == i), (parti == i))
         s <- sum(s[upper.tri(s)])
         s.min <- sum(sort(dis[upper.tri(dis)])[1:sum(parti == i)])
         s.max <- sum(sort(dis[upper.tri(dis)], decreasing = TRUE)[1:sum(parti == i)])          
         c.index <- c(c.index, (s - s.min)/(s.max - s.min))
      }
      q <- sum(c.index)
   } else if (method == "g2") {       
      q <- cluster.stats(dis, clustering = parti, G2 = TRUE)$g2  #large is good
   } else if (method == "g3") {       
      q <- cluster.stats(dis, clustering = parti, G3 = TRUE)$g3  #small is good
   } else if(method == "igp") {              #larger is good
      if (is.null(X))
         stop("data matrix is needed for 'igp' measure")
      if (any(is.na(X))) {
      	 X <- .myImpute(X, maxmiss)
      }	 
      Centroids <- do.call(cbind, lapply(unique(parti), function(i) rowMeans(X[, which(parti == i)])))
      colnames(Centroids) <- paste("centroid", unique(parti), sep="")
      if (is.null(rownames(X))) {
         rownames(X) <- paste("gene", 1:nrow(X), sep="")
         rownames(Centroids) <- paste("gene", 1:nrow(X), sep="")
      }
      if (is.null(colnames(X))) {
         rownames(X) <- paste("s", 1:ncol(X), sep="")
      }
      q <- mean(IGP.clusterRepro(X, Centroids)$IGP) 
   } else {
   	L <- cluster.stats(dis, clustering = parti, ...)
   	if (any(names(L) == method)) {
   	   q <- L[[which(names(L) == method)]]
   	} else {
   	   stop("given input 'method' does not exist")
   	}
   }	  	 	   	  
   return(q)

}
