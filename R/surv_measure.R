surv_measure <-
function(x, surv.time, status, method = "BIC")
{
   if (length(x) != length(surv.time) || length(status) != length(surv.time))  
      stop("input arguments dimension mismatch")
   if (!any(status %in% c(0, 1)))
      stop("Status vector must be zeros and ones")
   if (!is.vector(surv.time, mode = "numeric"))
      stop("'surv.time' should be numeric vector")   
   if (sum(method %in% c("AIC", "BIC")) == 0)
      stop("'method' must be 'AIC' or 'BIC'")
   data <- data.frame(st = surv.time, event = status, var = as.factor(x))
   fit <- coxph(Surv(st, event) ~ var, data)
   p <- length(unique(x))
   n <- length(x)
   d <- sum(status)  
   if (method == "AIC") {
      aic <- - 2*fit$loglik[2] +  2*(p+2)
      aic <- (aic + (2 * (p + 2) * (p + 3))/(n - p - 3))
      return(-aic)
   } else if (method == "BIC") {
      bic <- (- 2*fit$loglik[2] +  p * log(d))
      return(-bic)
   }
}
