cluster_pred <-
function(X, partition, surv.time = NULL, status = NULL, te.index, minclus = 4, 
                te.surv.time = NULL, te.status = NULL, method = "conc", maxmiss = 30, plot.it = FALSE, ...)           
{
      if (class(X) == "ExpressionSet" | class(X) == "eSet") {
         X <- exprs(X)
      } 
      if (length(partition) != ncol(as.matrix(X)[, -te.index]))
         stop("Input arguments dimension mismatch")
      if (!any(method %in% c("ward", "conc")))
         stop("'method' should be either 'ward' or 'conc'")   
      if (class(X) != "dist" & !isSymmetric(as.matrix(X))) {       
         if (any(is.na(X))) {
         	X <- .myImpute(X, maxmiss, ...)
         } 	
      }
      if (class(X) == "dist" | isSymmetric(as.matrix(X))) {
      	 if (any(is.na(as.matrix(X))))
      	    stop("Distance matrix must not has missing values")
      }  
      tr.index <- seq(ncol(X))[-te.index]
      X <- as.matrix(X)         
      if (method == "ward") {  ################################## Using the Ward distance approach
         best.index <- lapply(unique(partition), function(x) which(partition == x))
         qua <- matrix(NA, length(te.index), length(best.index))
         for(t in 1:length(te.index)) {
            d <- c()
            for(t1 in 1:length(best.index)) {         
               s.index <- c(te.index[t], tr.index[best.index[[t1]]])
               if (!isSymmetric(X)) {
                  D1 <- as.dist(1 - cor(X[, s.index]))
                  D2 <- as.dist(1 - cor(X[, s.index[-1]]))
               } else {
                  D1 <- as.dist(X[, s.index][s.index, ])
                  D2 <- as.dist(X[, s.index[-1]][s.index[-1], ])
               }
               qua[t, t1] <-  max(hclust(D1, method = "ward")$height) - max(hclust(D2, method = "ward")$height)
            }
         }
         pred <- apply(qua, 1, function(x) unique(partition)[which.min(x)])
         ############################tiny clusters pruning process
         chosen <- rep(1, length(pred))
         Do <- any(table(pred) < minclus) 
         while (Do) {
               freq <- as.vector(table(pred))
               g <- sort(unique(pred))[which.min(freq)]
               for(j in which(pred == g)) {
                  chosen[j] <- chosen[j] + 1
                  opt <- order(qua[j, ])[chosen[j]]
                  pred[j] <- unique(partition)[opt]
               }
               if (all(table(pred) >= minclus)) {
                  Do <- FALSE
               } 
         }      
      } else if (method == "conc") {      ################### Using Harrold's concordance approach
         if (is.null(surv.time) | is.null(status)) 
            stop("'surv.time' and 'status' must be supplied when method set to 'conc'")
         if ((length(partition) != length(tr.index)) | (length(surv.time) != length(status)))
            stop("Input arguments dimensions mismatch")
         if (!any(status %in% c(0, 1)))
            stop("'status' vector must be zeros and ones")
         if (!is.vector(surv.time, mode = "numeric"))
            stop("'surv.time' must be a numeric vector")     
         best.index <- lapply(unique(partition), function(x) which(partition == x))          
         K <- minclus
         if (!isSymmetric(X)) {
            ord <- lapply(te.index, function(x) order((cor(X[, c(x, tr.index)])[1, -1])[1:K], decreasing = TRUE))
         } else {
            ord <- lapply(te.index, function(x) order(X[x, tr.index])[1:K])
         }
         U <- unique(unlist(ord))
         val <- matrix(NA, length(U), length(best.index))
         for(u in 1:length(U)) {
            psuedo.test.st <- surv.time[U[u]]
            psuedo.test.event <- status[U[u]]
            for(t1 in 1:length(best.index)) {
               own.st <- surv.time[best.index[[t1]]]
               own.event <- status[best.index[[t1]]]
               ref.st <- surv.time[-best.index[[t1]]]
               ref.event <- status[-best.index[[t1]]]
               st.sim <- 0
               n <- 1 
               for(i in 1:length(own.st)) {
                  ref.st <- c(ref.st, own.st[-i])
                  ref.event <- c(ref.event, own.event[-i]) 
                  for(j in 1:length(ref.st)) {
                     pre.st.sim <- st.sim	
                     if (ref.event[j] == 1) {
                        if ((psuedo.test.event + own.event[i] == 0) & all(c(psuedo.test.st, own.st[i]) > ref.st[j])) { 
                           st.sim <- st.sim + 1
                        } else if (psuedo.test.event + own.event[i] == 2) {
                           st.sim <- st.sim + all(c(psuedo.test.st, own.st[i]) >= ref.st[j]) + 
                                         all(c(psuedo.test.st, own.st[i]) < ref.st[j])
                           st.sim <- st.sim - (psuedo.test.st >= ref.st[j] & own.st[i] < ref.st[j]) - 
                                         (own.st[i] >= ref.st[j] & psuedo.test.st < ref.st[j])
                        } else {
                           st.sim <- st.sim - psuedo.test.event*(ref.st[j] >= psuedo.test.st & ref.st[j] < own.st[i]) - 
                                         own.event[i]*(ref.st[j] >= own.st[i] & ref.st[j] < psuedo.test.st)
                        }
                     } else {
                        if ((psuedo.test.event + own.event[i] == 2) & all(ref.st[j] >= c(psuedo.test.st, own.st[i]))) {
                           st.sim <- st.sim + 1
                        }
                     }   
                     if (pre.st.sim != st.sim) {
                        n <- n + 1
                     }	  
                  }
               }
               val[u, t1] <- st.sim / n
            }
         }
         ord <- do.call(rbind, ord)
         qua <- array(NA, dim=c(length(te.index), length(best.index), K)) 
         for(t in 1:length(te.index)) {
            for(k in 1:K) {
               qua[t, , k] <- val[which(ord[t, k] == U), ] / k 
            }
         }
         pred <- unlist(lapply(1:length(te.index), function(x) unique(partition)[which.max(rowSums(qua[x, , ]))]))
         chosen <- rep(1, length(pred))
         ############################tiny clusters pruning process
         Do <- any(table(pred) < minclus) 
         while (Do) {
              freq <- as.vector(table(pred))
              g <- sort(unique(pred))[which.min(freq)]
              for(j in which(pred == g)) {
                 chosen[j] <- chosen[j] + 1
                 opt <- order(unlist(lapply(1:dim(qua)[2], function(i) sum(qua[j, i, ]))), 
                              decreasing = TRUE)[chosen[j]]
                 pred[j] <- unique(partition)[opt]
             }
             if (all(table(pred) >= minclus)) {
                Do <- FALSE
             } 
          }
      }
      if (!is.null(te.surv.time) & !is.null(te.status)) {
         St <- .KM_fit(te.surv.time, te.status, pred) 	
         pval <- NULL
         if (length(unique(pred)) >= 2) {
            pval <- pvalue(surv_test(Surv(te.surv.time, te.status) ~ as.factor(pred)))
         } 
         if (plot.it) {
            .plot_KM(St)         
         }
         Res <- list(St = St, pvalue = pval)
      } else {
         Res <- pred
      }       	   
    return(Res)   	  	  
}
