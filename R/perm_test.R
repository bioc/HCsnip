perm_test <-
function(partitions, surv.time, status, score1 = NULL, score2, method = "BIC", 
                      nperm = 1000) 
{

   if (!any(status %in% c(0, 1)))
      stop("Status vector must be zeros and ones")
   if (!is.vector(surv.time, mode = "numeric"))
      stop("'surv.time' must be a numeric vector")   
   if ((ncol(partitions) != length(surv.time)) || (ncol(partitions) != length(status)))
      stop("input arguments dimension mismatch") 
   if (is.null(score1)) {
   	  score1 <- apply(partitions, 1, function(x) surv_measure(x, surv.time, status, method))
   }	      
   ord <- order(score1, decreasing = TRUE)
   ranksum <- (1:length(ord)) + unlist(lapply(ord, function(x) which(x == order(score2, decreasing = TRUE))))
   best <- partitions[ord[which.min(ranksum)], ]
   obs.p <- pvalue(surv_test(Surv(surv.time, status) ~ as.factor(best)))
   perm.p <- c() 
   for (i in 1:nperm) {
       shaf <- sample(1:length(status), length(status), replace = FALSE)
       st <- surv.time[shaf]
       event <- status[shaf]
       L <- apply(partitions, 1, function(x) surv_measure(x, st, event, method))
       ord <- order(L, decreasing = TRUE) 
       ranksum <- (1:length(ord)) + unlist(lapply(ord, function(x) which(x == order(score2, decreasing = TRUE))))
       best.p <- partitions[ord[which.min(ranksum)], ]
       perm.p <- c(perm.p, pvalue(surv_test(Surv(st, event) ~ as.factor(best.p))))
   }
   cat("Corrected p-value of the best partition is",  sum(perm.p <= obs.p)/nperm, "\n")

 return(list(obs.p = obs.p, perm.p = perm.p, best = best))
   
}
