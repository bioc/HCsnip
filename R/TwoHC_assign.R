TwoHC_assign <-
function(X, index1, index2, new.X, dis.method = "cor", link.method = "ward", minclus = 4,
             maxmiss = 30, surv.time, status, method1 = "BIC", method2 = "g2")
{
   if (class(X) == "ExpressionSet" | class(X) == "eSet") {
      X <- exprs(X)
   }
   if ((length(surv.time) != ncol(X)) || (ncol(X) != length(c(index1, index2)))) 
      stop("Input arguments dimension mismatch")
   #####group 1
   res1 <- HCsnipper(X = X[, index1], dis.method = dis.method, link.method = link.method, minclus = minclus, maxmiss = maxmiss)
   dis.1 <- res1$dis
   H.1 <- res1$hc
   cl.1 <- res1$partitions
   L1.1 <- apply(cl.1, 1, function(x) surv_measure(x, surv.time[index1] , status[index1], method1))
   L2.1 <- apply(cl.1, 1, function(x) measure(x = x, dis = dis.1, type = method2))
   ord <- order(L1.1, decreasing = TRUE)
   ranksum <- (1:length(ord)) + unlist(lapply(ord, function(x) which(x == order(L2.1, decreasing = TRUE))))	     
   best1 <- cl.1[ord[which.min(ranksum)], ]
   #####group 2
   res2 <- HCsnipper(X = X[, index2], dis.method = dis.method, link.method = link.method, minclus = minclus, maxmiss = maxmiss)
   dis.2 <- res2$dis
   H.2 <- res2$hc
   cl.2 <- res2$partitions
   L1.2 <- apply(cl.2, 1, function(x) surv_measure(x, surv.time[index2] , status[index2], method1))
   L2.2 <- apply(cl.2, 1, function(x) measure(x = x, dis = dis.2, type = method2))
   ord <- order(L1.2, decreasing = TRUE)
   ranksum <- (1:length(ord)) + unlist(lapply(ord, function(x) which(x == order(L2.2, decreasing = TRUE))))	      
   best2 <- cl.2[ord[which.min(ranksum)], ]
   ##########################################
   if (class(X) == "ExpressionSet" | class(X) == "eSet") {
      new.X <- exprs(new.X)
   }
   if (any(is.na(new.X))) {
      new.X <- .myImpute(new.X, maxmiss)
   }   
   index <- c(list(index1), list(index2))
   best <- c(list(best1), list(best2))
   res <- t(apply(new.X, 2, function(x) .test_pred(best, x, surv.time, status, X, index, minclus)))
   colnames(res) <- c("which.HC", "which.cluster.hc1", "which.cluster.hc2")
   m1 <- data.frame(L1.1, L2.1)
   m2 <- data.frame(L1.2, L2.2)
   colnames(m1) <- colnames(m2) <- c(method1, method2)
   result <- list(hc1 = H.1, hc2 = H.2, partitions.hc1 = cl.1, partitions.hc2 = cl.2, best.hc1 = best1, best.hc2 = best2, 
                  score.hc1 = m1, score.hc2 = m2, Assign = res, surv.time = surv.time, status = status,
                  index1 = index1, index2 = index2, new.X = new.X, X = X, method1 = method1, method2 = method2,
                  minclus = minclus)
   
 return(result)
         
}
