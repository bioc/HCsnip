EnvioPlot <- function(dat, method = "knn", x,
               horizontal = FALSE, col = NULL, names = NULL, ...)
{
    if (class(dat) == "ExpressionSet" | class(dat) == "eSet") {
       dat <- exprs(dat)
    }
    if (any(is.na(dat))) {
       dat <- .myImpute(dat, maxmiss = 30)
    }
    n <- length(unique(x))
    at <- 1:n
    if (!is.null(col)) {
       if (length(col) != n)
          stop("The number of color to use not equal the number clusters")
    } else {
       col <- 2:(n+1)
    }
    H <- hdEntropy(t(dat), method = method, ...)
    upper <- vector(mode = "numeric", length = n)
    lower <- vector(mode = "numeric", length = n)
    q1 <- vector(mode = "numeric", length = n)
    q3 <- vector(mode = "numeric", length = n)
    med <- vector(mode = "numeric", length = n)
    base <- vector(mode = "list", length = n)
    height <- vector(mode = "list", length = n)
    baserange <- c(Inf, -Inf)
    args <- list(display = "none")
    datas <- lapply(sort(unique(x)), function(i) H[x == i])
    for (i in 1:n) {
        data <- datas[[i]]
        data.min <- min(data)
        data.max <- max(data)
        q1[i] <- quantile(data, 0.25)
        q3[i] <- quantile(data, 0.75)
        med[i] <- median(data)
        iqd <- q3[i] - q1[i]
        upper[i] <- min(q3[i] + 1.5 * iqd, data.max)
        lower[i] <- max(q1[i] - 1.5 * iqd, data.min)
        est.xlim <- c(min(lower[i], data.min), max(upper[i], 
            data.max))
        smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
            args))
        hscale <- 0.4/max(smout$estimate) 
        base[[i]] <- smout$eval.points
        height[[i]] <- smout$estimate * hscale
        t <- range(base[[i]])
        baserange[1] <- min(baserange[1], t[1])
        baserange[2] <- max(baserange[2], t[2])
    }
    ylim <- baserange
    if (n == 1) { 
       xlim <- at + c(-0.5, 0.5)
    } else { 
       xlim <- range(at) + min(diff(at))/2 * c(-1, 1)
    }
    if (is.null(names)) {
        label <- paste("cluster", sort(unique(x)), sep="")
    } else {
        label <- names
    }
    boxwidth <- 0.05 
    plot.new()
    if (!horizontal) {
       plot.window(xlim = xlim, ylim = ylim)
       axis(2)
       axis(1, at = at, labels = label)
       box()
       for (i in 1:n) {
           polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
                c(base[[i]], rev(base[[i]])), col = col[i], border = "black", 
                lty = 1, lwd = 1)
           lines(at[c(i, i)], c(lower[i], upper[i]), lwd = 1, lty = 1)
           rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, q3[i], col = "black")
           points(at[i], med[i], pch = 19, col = "white")
       }
    } else {
        plot.window(xlim = ylim, ylim = xlim)
        axis(1)
        axis(2, at = at, labels = label)
        box()
        for (i in 1:n) {
            polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
                rev(at[i] + height[[i]])), col = col[i], border = "black", 
                lty = 1, lwd = 1)
            lines(c(lower[i], upper[i]), at[c(i, i)], lwd = 1, lty = 1)
            rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + boxwidth/2, col = "black")
            points(med[i], at[i], pch = 19, col = "white")
        }
    }
    invisible(list(upper = upper, lower = lower, median = med, 
        q1 = q1, q3 = q3))
    title(main = "Cluster Entropies")
    title(ylab = "Entropy")
    return(H)
}








