.test_pred <- function(best, t.em, surv.time, status, em, index, k)
{
  
     if (!any(status %in% c(0, 1)))
        stop("Status vector must be zeros and ones")
     if (!is.vector(surv.time, mode = "numeric"))
        stop("'surv.time' should be numeric vector")
     cl.time <- cl.event <- c()   
     cl.index <- rep(NA, length(best))
     for(p in 1:length(best)) {
        tr.index <- index[[p]]
        st <- surv.time[tr.index]
        event <- status[tr.index]
        best.index <- lapply(unique(best[[p]]), function(x) which(best[[p]] == x))
        U <- order(cor(cbind(t.em, em[, tr.index]))[1, -1], decreasing = TRUE)[1:k]
        val <- matrix(NA, length(U), length(best.index))
        for(u in 1:length(U)) {
           psuedo.test.st <- st[U[u]]
           psuedo.test.event <- event[U[u]]
           for(t1 in 1:length(best.index)) {
              own.st <- st[best.index[[t1]]]
              own.event <- event[best.index[[t1]]]
              ref.st <- st[-best.index[[t1]]]
              ref.event <- event[-best.index[[t1]]]
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
        val <- val / length(best.index)  ##### take into account the number 
                                         ##### of clusters in the best partitions of treaments
        #pred[p] <- max(apply(val, 2, function(x) as.vector(x %*% (1/(1:k)))))
        cl.index[p] <- which.max(apply(val, 2, function(x) as.vector(x %*% (1/(1:k)))))
        cl.time <- c(cl.time, st[best[[p]] == cl.index[p]]) 
        cl.event <- c(cl.event, event[best[[p]] == cl.index[p]])
     }   
     cls <- c(rep(1, sum(best[[1]] == cl.index[1])), rep(0, sum(best[[2]] == cl.index[2])))
     fit <- coxph(Surv(cl.time, cl.event) ~ factor(cls))
     if (coef(fit) < 0) { 
        result <- c(1, cl.index)
     } else {
     	result <- c(2, cl.index)
     }	      	       	 
  return(result)
}


#################################################################
### function to compute Kaplan-Meier estimates
###
.KM_fit <- function(times, censor, group)
{
     cens <- gp <- tt <- NULL
     j <- 0
	 for(i in times) {
	     j <- j + 1
	     tt <- c(tt, i)
	     cens <- c(cens, rep(1, length(i)))
	     cens[length(cens)] <- censor[j]
	     gp <- c(gp, rep(group[j], length(i)))
	 }
	 times <- tt
	 censor <- cens
	 group <- gp
	 rm(tt, cens, gp)
     group <- as.numeric(group)
     if (!is.vector(times, mode = "numeric"))
        stop("'surv.time' must be a numeric vector")
     if (any(times < 0))
        stop("negative times")
     if (length(group) == 1)
        group <- rep(1, length(times))
     if (length(censor) == 1)
        censor <- rep(1, length(times))
     if (length(group) != length(times) || length(times) != length(censor))
	    stop("All vectors must be the same length")
     if (!any(censor %in% c(0, 1)))
	    stop("'status' vector must be zeros and ones")
     
     # calculate unique times of events
     #
     ss <- cens <- gg <- ff <- cg <- NULL
     for(i in unique(group)) {
	    n <- sum(group == i)
	    o <- order(times[group == i])
	    cens2 <- censor[group == i][o]
	    times2 <- times[group == i][o]
	    f2 <- as.vector(table(list(1 - cens2, times2)))
	    if (!any(cens2 == 0) || !any(cens2 == 1)) {
		   ff2 <- rep(0, 2 * length(f2))
		   ff2[seq(1, 2 * length(f2) - 1, by = 2)] <- f2
		   f2 <- ff2
		}
	    l2 <- length(f2)
	    g2 <- rep(i,l2)
	    c2 <- rep(c(1,0),l2/2)
	    s2 <- vector(mode = "numeric", length = l2)
	    s2[seq(1, l2 - 1, by = 2)] <- unique(times2)
	    s2[seq(2, l2, by = 2)] <- unique(times2)
	    ss <- c(ss, s2)
	    cens <- c(cens, c2)
	    gg <- c(gg, g2)
	    ff <- c(ff, f2)
	    cg <- c(cg, f2[l2] == 0)
	 }
     # calculate survivor curve
     #
     j <- (ff != 0) | (cens == 1)
     k <- seq(1:length(ff))
     n <- sum(j)
     t <- cc <- f <- g <- tt <- vector(mode = "numeric", length = n)
     t <- ss[k * j]
     cc <- cens[k * j]
     f <- ff[k * j]
     g <- gg[k * j]
     v <- r <- s <- rep(0, n)
     for(i in unique(group)) {
	    tt[n:1] <- cumsum(f[n:1] * (g[n:1] == i))
	    r <- r + tt * (g == i)
	    tmp <- cumsum(log(ifelse(tt == f | tt == 0, 1, (tt - f) / tt)) * (g == i) * cc)
	    s <- s + exp(ifelse(is.na(tmp), 0, tmp)) * (g == i) * (tt != f)
	    tmp <- cumsum(f / tt / (tt - f) * (g == i) * cc)
	    v <- v + ifelse(is.na(tmp), 0, tmp) * (g == i)
	 }
     v <- s * s * v 
     # set up matrix to return
     #
     m <- NULL
     for(i in unique(g)) m <- c(m, sum(g == i))
     j <- seq(1, length(g))
     j[cumsum(m)] <- cg * j[cumsum(m)]
     z <- cbind(t[j], g[j], r[j], s[j], v[j])
     dp <- rep(TRUE, dim(z)[1])
     for(i in 2:dim(z)[1]) if (all(z[i,] == z[i-1,], na.rm = TRUE)) dp[i] <- FALSE
     z <- z[dp, ]
     colnames(z) <- c("Time", "Cluster", "At risk", "S(t)", "Var(S)") 
     rownames(z) <- paste(rep("",dim(z)[1]))
   return(z)  
}


### function to plot survivor curves output from km
###
.plot_KM <- function(z)
{
     surv <- z[, 4]
     times <- z[, 1]
     group <- z[, 2]
     kms <- ttt <- NULL
     add <- FALSE
     # plot one curve per group
     #
     k <- ln <- lt <- 0
     col <- 1:length(unique(group))
     pch <- c()
     nn <- 1  
     for(i in unique(group)) {
	    lt <- lt %% 4+1
          pch <- c(pch, lt)
	    if (length(group) > 1 && length(unique(group)) > 1) {
		   j <- (cumsum(group == i) + ln) * (group == i)
		   s <- surv[j]
		   if (!missing(times)) t <- times[j]
		   ln <- ln + sum(group == i)
	    } else {
		   s <- surv
		   t <- times
	    }
	    n <- 2 * length(s) - 1
	    # prepare survivor function
	    km <- rep(0, n)
	    km[seq(1, n, by = 2)] <- s
	    km[seq(2, n - 1, by = 2)] <- s[1:(length(s) - 1)]
	    km <- c(1, 1, km)
	    # prepare times
     	tt <- rep(0,n)
    	tt[seq(1, n, by = 2)] <- t
	    tt[seq(2, n - 1, by = 2)] <- t[1:(length(s) - 1)]
	    tt <- c(0, tt, tt[length(tt)])
	    # plot
	    if (add)  
	       lines(tt, km, lty = lt, pch = lt, col = col[nn])
	    else {
		   xlim <- c(min(times) - 1, max(times + 1))
		   plot(tt, km, type = "s", xlim = xlim, ylim = c(0, 1), ylab = "Survival probability", xlab = "Time to event",
		        main = "Kaplan-Meier curve", lty = lt, pch = lt, col = col[nn])
	    }
	    add <- TRUE
	    ttt <- c(ttt, tt)
	    kms <- c(kms, km)
	    nn <- nn + 1
     }
     legend("topright", cex = 1, pch = pch[order(unique(group))], col = col[order(unique(group))], 
                           legend = paste("cluster", sort(unique(group)), sep = ""), ncol = 1)
     invisible(cbind(ttt, kms))
}


##################### function to impute missing values
##
.myImpute <- function(X, maxmiss = 30, ...)
{
	
    countMissing <- function(x, maxMissing) 
    {
          if (length(x[is.na(x)]) <= maxMissing) 
               return(TRUE)
            return(FALSE)
    }
    nummis <- length(X[is.na(X)])
    if (nummis > 0) {
       allowed <- ncol(X) * maxmiss/100
       filter <- apply(X, 1, countMissing, allowed)
       X <- X[filter, ]
    }
    nummis <- length(X[is.na(X)])
    if (nummis > 0) {
       if (!exists("k")) 
          k <- 10
       new.k <- ceiling((1 - maxmiss/100) * ncol(X))
       new.k <- min(new.k, k)
       if (new.k != 10) 
          cat("Changing impute.knn parameter k from", k, "to", 
                  new.k, "due to small sample size.\n")
       if (new.k > 1) {
          nrows <- nrow(X)
          nrpi <- floor(nrows/100)
          resultsall <- c()
          for (i in 1:(100 - 1)) {
              datforimp <- X[((i - 1) * nrpi + 1):(nrpi * i), ]
              result <- impute::impute.knn(datforimp, k = new.k, ...)
              resultsall <- rbind(resultsall, result$data)
          }
          datforimp <- X[((100 - 1) * nrpi + 1):nrows, ]
          result <- impute::impute.knn(datforimp, k = new.k, ...)
          X <- rbind(resultsall, result$data)
       }
    }
    return(X)
}

