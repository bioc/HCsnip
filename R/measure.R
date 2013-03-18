measure <-
function(x, dis, dat = NULL, method = "g2", maxmiss = 30, ...) 
{
   if (length(x) != nrow(as.matrix(dis))) 
      stop("input arguments dimension mismatch")
   if (!isSymmetric(as.matrix(dis))) 
      stop("'dis' should be a square distance matrix or 'dist' class object")
   if (method == "c.index") {  #small is good 
      c.index <- c()
      dis <- as.matrix(dis)
      for(i in unique(x)) {
         s <- subset(dis, (x == i), (x == i))
         s <- sum(s[upper.tri(s)])
         s.min <- sum(sort(dis[upper.tri(dis)])[1:sum(x == i)])
         s.max <- sum(sort(dis[upper.tri(dis)], decreasing = TRUE)[1:sum(x == i)])          
         c.index <- c(c.index, (s - s.min)/(s.max - s.min))
      }
      q <- sum(c.index)
   } else if (method == "g2") {       
      q <- cluster.stats(dis, clustering = x, G2 = TRUE)$g2  #large is good
   } else if (method == "g3") {       
      q <- cluster.stats(dis, clustering = x, G3 = TRUE)$g3  #small is good
   } else if(method == "IGP") {              #larger is good
      if (is.null(dat))
         stop("data matrix is needed for 'IGP' measure")
      if (any(is.na(dat))) {
      	 dat <- .myImpute(dat, maxmiss)
      }	 
      Centroids <- do.call(cbind, lapply(unique(x), function(i) rowMeans(dat[, which(x == i)])))
      colnames(Centroids) <- paste("centroid", unique(x), sep="")
      if (is.null(rownames(dat))) {
         rownames(dat) <- paste("gene", 1:nrow(dat), sep="")
         rownames(Centroids) <- paste("gene", 1:nrow(dat), sep="")
      }
      if (is.null(colnames(dat))) {
         rownames(dat) <- paste("s", 1:ncol(dat), sep="")
      }
      q <- mean(IGP.clusterRepro(dat, Centroids)$IGP) 
   } else {
   	L <- cluster.stats(dis, clustering = x, ...)
   	if (any(names(L) == method)) {
   	   q <- L[[which(names(L) == method)]]
   	} else {
   	   stop("given input 'method' does not exist")
   	}
   }	  	 	   	  
   return(q)

}
