TwoHC_perm <-
function(TwoHC, nperm = 1000)
{

   wrapper <- function(surv.time, status, index, L2.1, L2.2, cl.1, cl.2, X, 
                       new.X, type1, minclus)
   {
       shaf.1 <- sample(1:length(index[[1]]), length(index[[1]]), replace = FALSE)
       shaf.2 <- sample(1:length(index[[2]]), length(index[[2]]), replace = FALSE)
       #####group 1
       L1.1 <- apply(cl.1, 1, function(x) surv_measure(x, surv.time[index[[1]]][shaf.1], 
                           status[index[[1]]][shaf.1], type1))                     
       ord <- order(L1.1, decreasing = TRUE)
       ranksum <- (1:length(ord)) + unlist(lapply(ord, function(x) which(x == order(L2.1, decreasing = TRUE))))
       best.1 <- cl.1[ord[which.min(ranksum)], ]                
       #####group 2
       L1.2 <- apply(cl.2, 1, function(x) surv_measure(x, surv.time[index[[2]]][shaf.2], 
                          status[index[[2]]][shaf.2], type1))
       ord <- order(L1.2, decreasing = TRUE)
       ranksum <- (1:length(ord)) + unlist(lapply(ord, function(x) which(x == order(L2.2, decreasing = TRUE))))
       best.2 <- cl.2[ord[which.min(ranksum)], ]                   
       ##########################################
       best <- c(list(best.1), list(best.2))
       res <- .test_pred(best, new.X, surv.time, status, X, index, minclus)
       best.index.1 <- which(best.1 == res[2])
       best.index.2 <- which(best.2 == res[3])
       cls <- rep(0, length(best.index.1) + length(best.index.2))
       cls[1:length(best.index.1)] <- 1
       sel <- c((index[[1]][shaf.1])[best.index.1], (index[[2]][shaf.2])[best.index.2]) 
       fit <- coxph(Surv(surv.time[sel], status[sel]) ~ cls)
    return(coef(fit))
   }
   ######
   ##
   grp <- TwoHC$Assign
   new.X <- TwoHC$new.X
   obs.betas <- c()
   perm.betas <- matrix(NA, nperm, ncol(new.X)) 
   for(i in 1:ncol(new.X)) {
      best.index.1 <- which(TwoHC$best.hc1 == grp[i, 2])
      best.index.2 <- which(TwoHC$best.hc2 == grp[i, 3])      
      cls <- rep(0, length(best.index.1) + length(best.index.2))
      cls[1:length(best.index.1)] <- 1
      sel <- c(TwoHC$index1[best.index.1], TwoHC$index2[best.index.2]) 
      fit <- coxph(Surv(TwoHC$surv.time[sel], TwoHC$status[sel]) ~ cls)
      obs.betas <- c(obs.betas, coef(fit))
      index <- c(list(TwoHC$index1), list(TwoHC$index2))
      surv.time <- TwoHC$surv.time
      status <- TwoHC$status
      X <- TwoHC$X
      L2.1 <- TwoHC$score.hc1[, 2]
      L2.2 <- TwoHC$score.hc2[, 2]
      cl.1 <- TwoHC$partitions.hc1
      cl.2 <- TwoHC$partitions.hc2
      type1 <- TwoHC$method1
      minclus <- TwoHC$minclus   
      for(b in 1:nperm) {
         perm.betas[b, i] <- wrapper(surv.time, status, index, L2.1, L2.2, cl.1, cl.2, X, 
                           new.X[, i], type1, minclus)
      }                     
   }
   r <- unlist(lapply(1:ncol(new.X), function(i)  
             sum(abs(obs.betas[i]) > abs(perm.betas[which(!is.na(perm.betas[, i])), i])) / sum(!is.na(perm.betas[, i]))))
   rr <- unlist(lapply(1:ncol(new.X), function(i)  
             exp(-abs(obs.betas[i])) / exp(-median(abs(perm.betas[which(!is.na(perm.betas[, i])), i])))))
   ##
   tooSmall <- function(x) {      ##used to avoid to small number
   	    for(i in 1:length(x)) {
   	       x[i] <- max(0.015, x[i])
   	    }
   	   return(x) 
   }	       	          
   gm.obs <- exp(log(mean(tooSmall(rr))))          
   z <-  unlist(lapply(1:ncol(new.X), function(i) exp(-median(abs(perm.betas[which(!is.na(perm.betas[, i])), i])))))
   gm.null <- rep(NA, nperm)
   for(i in 1:nperm) {
   	  a <- exp(-abs(perm.betas[i, ]))
   	  gm.null[i] <- exp(log(mean(tooSmall(a[!is.na(a)] / z[!is.na(a)]))))
   }
   pval <- sum(gm.obs >= gm.null) / nperm	            
   result <- list(Obs.betas = obs.betas, Perm.betas = perm.betas, Ranks = r, RiskRatios = rr, Pvalue = pval)
   
  return(result)
}










