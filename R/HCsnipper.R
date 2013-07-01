HCsnipper <-
function(X, hc = NULL, dis = NULL, dis.method = "cor", link.method = "ward", minclus = 4, maxmiss = 30, ...) 
{
  	
  D <- NULL
  if (is.null(hc)) {
     if (is.null(dis)) {
        if (any(is.na(X))) {
           X <- .myImpute(X, maxmiss, ...)
        }   	
        if (dis.method == "cor") {
           hc <- hclust(as.dist(1 - cor(X)), method = link.method)
           D <- 1 - cor(X)
        } else {
           hc <- hclust(dist(t(X), method = dis.method), method = link.method)
           D <- as.matrix(dist(t(X)))
        }
     } else {
        if (class(dis) == "dist") {
           hc <- hclust(dis, method = link.method)
           D <- as.matrix(dis)
        } else {
           hc <- hclust(as.dist(dis), method = link.method)
           D <- dis
        }
     } 
  }
  if (!exists("X", mode = "numeric"))
     X <- NULL
  if (class(hc) != "hclust") 
     stop("'hc' should be hclust class object") 
  CrossOver <- function(CLUSTER, Mother, Child) 
   {
      Size <- sapply(CLUSTER, function(x) nrow(x))
      Cluster <- matrix(0 , nrow = prod(Size), ncol = length(Mother))
      for(t in 1:max(Mother)) {
         cl <- CLUSTER[[t]]
         if (t < max(Mother)) {  
            Step <- seq(1, prod(Size), prod(Size[-(1:t)]))
         }else {
            Step <- seq(prod(Size))
         }  
         if (t == 1) {
            for(r in 1:length(Step)) {
               Cluster[Step[r]:(Step[r] + prod(Size[-c(1:t)])-1), Child[[t]]] <- 
                   do.call(rbind, rep(list(cl[r, Child[[t]]]), prod(Size[-(1:t)])))
            }
         } else {
            for(i in 1:prod(Size[1:(t-1)])) {
               for(r in 1:Size[t]) {
                  Add <- unique(Cluster[Step[r]:(Step[r] + prod(Size[-c(1:t)])-1), -Child[[t]]])
                  Stack <- as.vector(do.call(rbind, rep(list(cl[r, Child[[t]]]), prod(Size[-(1:t)]))))
                  U <- Stack
                  n <- 1
                  for(u in unique(U)) {
                     Stack[which(U == u)] <-  max(Add) + n
                     n <- n+1 
                  } 
                  Cluster[Step[r]:(Step[r] + prod(Size[-c(1:t)])-1), Child[[t]]] <- 
                      matrix(Stack, ncol=length(Child[[t]]))
               }
               Step <- Step[-c(1:Size[t])]
            }
         }
      }
      Cluster <- unique(Cluster)
      return(Cluster)
   }
   Clus.index <- vector("list", nrow(hc$merge))
   for(t in 1:(nrow(hc$merge))) {
      if (all(hc$merge[t, ]<0)) {
         Clus.index[t] <- list(abs(hc$merge[t, ]))
      } else if (any(hc$merge[t,]<0)) {
         Clus.index[t] <- list(c(abs(hc$merge[t, which.min(hc$merge[t, ])]), 
              unlist(Clus.index[hc$merge[t, which.max(hc$merge[t, ])]])))
      } else {
         Clus.index[t] <- list(c(unlist(Clus.index[hc$merge[t, 1]]), 
              unlist(Clus.index[hc$merge[t, 2]])))
      }
   }
   Merge <- data.frame(hc$merge, hc$height)
   L <- mean(tail(Merge[, 3], 2))
   Mother <- as.vector(cutree(hc, h = L))
   Child <- lapply(1:max(Mother), function(x) which(Mother == x))
   Child <- lapply(1:max(Mother), function(x) hc$order[hc$order %in% Child[[x]]])
   p <- sapply(Clus.index, function(x) sapply(Child, function(y) all(x %in% y)))
   CUT <- list()
   for(i in 1:max(Mother)) {
      CUT <- c(CUT, list(which(p[i, ] == TRUE)))
   }
   cluster <- list()
   for(m in 1:max(Mother)) {   
      Ind <- CUT[[m]]
      W <- Clus.index[Ind]
      W <- W[order(sapply(W, function(x) length(x)))]
      C <- c()
      go <- TRUE
      j <- 1
      if (length(Child[[m]]) > 2) {
         while(j <= (length(Ind) - 1)) {
         	  L <- mean(Merge[rev(Ind)[j:(j + 1)], 3])
              cut <- as.vector(cutree(hc, h = L))
              C <- rbind(C, rep(0, length(Mother)))
              C[nrow(C), Child[[m]]] <- cut[Child[[m]]]
              j <- j + 1
         }
         cl <- c()
         for(j in 1:nrow(C)) {
         	if (any(table(C[j, ]) == 1)) {
         	   singl <- as.numeric(names(which(table(C[j, ]) == 1)))
         	   for(jj in 1:length(singl)) {
         	      idx <- which(C[j, ] == singl[jj])	
                  id <- which(sapply(W, function(x) (!any(idx %in% x) & !all(Child[[m]] %in% c(idx, x)))))
                  for(jjj in id) {  
                     if (length(unique(C[j, W[[jjj]]])) == 1) {
                        cl <- rbind(cl, C[j, ])
               	        cl[nrow(cl), idx] <- C[j, W[[jjj]]][1]
                     }
                  }
               }
            }
         }
         C <- rbind(C, cl)        	  
         C <- unique(C)
         for(i in 1:length(W)) {
            cl <- c()
            for(j in 1:(length(Ind) - 1)) {
               L <- mean(Merge[rev(Ind)[j:(j + 1)], 3])	
               cl <- rbind(cl , rep(0, length(Mother)))
               cut <- as.vector(cutree(hc, h = L))
               cl[nrow(cl), ] <- cut
               cl[nrow(cl), W[[i]]] <- max(cl[nrow(cl), -(W[[i]])], 0) + 1
            }
            for(y in 1:nrow(cl)) {
               U <- unique(cl[y, Child[[m]]])
               b <- cl[y, Child[[m]]]
               for(u in 1:length(U)) {
                  cl[y, Child[[m]][which(b == U[u])]] <- u
               }
            }
            cl[,  Child[[seq(max(Mother))[-m]]]] <- 0
            cl <- unique(cl)
            C <- rbind(C, cl)
         }
         C <- unique(C)
         no <- c()
         for(i in 1:nrow(C)) {
   	      a <- sapply(W[which(sapply(W, length) == 2)], function(x) length(unique(C[i, x])))
            if (any(a != 1)) {
               no <- c(no, i)
            }
         }
         C <- C[-no, ] 
         if (length(nrow(C)) != 0) {
            cluster <- c(cluster, list(unique(C)))
         } else {
            cluster <- c(cluster, list(C))
         }
      } else {
         a <- rep(0, length(Mother))
         a[Child[[m]]] <- 1
         a <- rbind(a, a)
         cluster <- c(cluster, list(a))
      }
   }   
   cluster <- lapply(1:max(Mother), function(x) rbind(cluster[[x]], Mother))
   cluster <- CrossOver(cluster, Mother, Child)
   for(y in 1:nrow(cluster)) {
      U <- unique(cluster[y, ])
      b <- cluster[y, ]
      for(u in 1:length(U)) {
         cluster[y, which(b == U[u])] <- u
      }
   }   
   cluster <- unique(cluster)  	 
   idx <- which(apply(cluster, 1, function(x) all(table(x) >= minclus)))
   cat("HC snipping is finished!", "\n")
   cat(paste(nrow(cluster), " unique partitions are found", sep=""), "\n")
   result <- list(partitions = cluster, hc = hc, dat = X, minclus = minclus, dis = D, 
   	                id = idx, dis.m = dis.method, link.m = link.method)
   	     
 
 return(result)
}
