RSF_eval <-
function(partition, surv.time, status, te.partition, te.surv.time, te.status, ...) 
{

   if (!any(status %in% c(0, 1)) || !any(te.status %in% c(0, 1)))
      stop("Status vector must be zeros and ones")
   if (!is.vector(surv.time, mode = "numeric") || !is.vector(te.surv.time, mode = "numeric"))
      stop("'surv.time' must be a numeric vector")     
   if (length(partition) != length(surv.time) || length(status) != length(surv.time))
      stop("input arguments dimensions mismatch")
   if (length(te.partition) != length(te.surv.time) || length(te.status) != length(te.surv.time))
      stop("input arguments dimensions mismatch")   
   fea.tr  <- matrix(0, length(partition), length(unique(partition)))
   for(j in 1:length(partition)) {
      fea.tr[j, partition[j]] <- 1
   }
   fea.te  <- matrix(0, length(te.partition), length(unique(partition)))
   for(j in 1:length(te.partition)) {
      fea.te[j, te.partition[j]] <- 1
   }
   if (any(surv.time == 0)) {
   	  surv.time <- surv.time + 1
   }
   if (any(te.surv.time == 0)) {
   	  te.surv.time <- te.surv.time + 1
   }
   colnames(fea.tr) <- colnames(fea.te) <- paste("cluster", 1:ncol(fea.tr), sep = "")
   tr <- data.frame(st = surv.time, event = status, fea.tr)
   te <- data.frame(st = te.surv.time, event = te.status, fea.te)   
   rf <- rfsrc(Surv(st, event)~., data = tr, ...)

 return(predict.rfsrc(rf, te)$err.rate)

}
