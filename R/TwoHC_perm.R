TwoHC_perm <-
function(TwoHC, nperm = 1000)
{

   wrapper <- function(q)
   {
       shaf <- sample(1:length(st), length(st), replace = FALSE)
       fit <- coxph(Surv(st[shaf], event[shaf]) ~ cls)
    return(coef(fit))
   }
   ###
   tooSmall <- function(x) {      ##used to avoid to small number
   	    for(i in 1:length(x)) {
   	       x[i] <- max(0.25, x[i])
   	    }
   	   return(x) 
   }
   ###
   ##
   tooBig <- function(x) {      ##used to avoid to small number
   	    for(i in 1:length(x)) {
   	       x[i] <- min(3, x[i])
   	    }
   	   return(x) 
   }
   ######
   ##
   grp <- TwoHC$Assign
   new.X <- TwoHC$new.X
   obs.betas <- rep(NA, ncol(new.X))
   perm.betas <- matrix(NA, nperm, ncol(new.X)) 
   for(i in 1:ncol(new.X)) {
      best.index.1 <- which(TwoHC$best.hc1 == grp[i, 2])
      best.index.2 <- which(TwoHC$best.hc2 == grp[i, 3])      
      cls <- rep(0, length(best.index.1) + length(best.index.2))
      cls[1:length(best.index.1)] <- 1
      sel <- c(TwoHC$index1[best.index.1], TwoHC$index2[best.index.2]) 
      obs.betas[i] <- coxph(Surv(TwoHC$surv.time[sel], TwoHC$status[sel]) ~ cls)$coef
      st <- TwoHC$surv.time[sel]          ##passed to the permutation
      event <- TwoHC$status[sel]
      perm.betas[, i] <- unlist(lapply(1:nperm, wrapper))
   }
   r <- unlist(lapply(1:ncol(new.X), function(i)  
             sum(abs(obs.betas[i]) > abs(perm.betas[which(!is.na(perm.betas[, i])), i])) / sum(!is.na(perm.betas[, i]))))
   rr <- unlist(lapply(1:ncol(new.X), function(i)  
             exp(-abs(obs.betas[i])) / exp(-median(abs(perm.betas[which(!is.na(perm.betas[, i])), i])))))
   rr <- tooBig(tooSmall(rr))
   obs.T <- mean(log(rr)) / sd(log(rr))    #observed test statistics
   z <-  unlist(lapply(1:ncol(new.X), function(i) exp(-median(abs(perm.betas[which(!is.na(perm.betas[, i])), i])))))
   null.T <- rep(NA, nperm)  ##test statistics from the null distribution
   for(i in 1:nperm) {
   	  a <- exp(-abs(perm.betas[i, ]))
   	  a <- tooBig(tooSmall(a[!is.na(a)] / z[!is.na(a)])) 
        null.T[i] <- mean(log(a)) / sd(log(a)) 
   }
   pval <- sum(obs.T >= null.T) / nperm	            
   result <- list(Obs.betas = obs.betas, Perm.betas = perm.betas, Ranks = r, RiskRatios = rr, Pvalue = pval)
   
  return(result)
}










