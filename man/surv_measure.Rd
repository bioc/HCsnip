\name{surv_measure}
\alias{surv_measure}

\title{
Cluster quality evaluation using follow-up data
}

\description{
Function to evaluate the overall quality of a given partition using follow-up data. A partition (clustering) is composed of non-overlapping clusters.
}

\usage{
surv_measure(parti, surv.time, status, method = "BIC")
}

\arguments{
  \item{parti}{A partition to be evaluated. }
  
  \item{surv.time}{A numeric vector contains follow-up time of patient's in \emph{x} }
  
  \item{status}{A binary vector contains survival status of patients in \emph{x}, normally 0=alive, 1=dead.}
  
  \item{method}{Type of partition evaluation measures to use for assessing the relationship between follow-up and a partition. Default is \emph{BIC}.}
  
}

\details{
This function fits a Cox model using follow-up data as response and cluster labels in the partition as covariate.  The likelihood from the fitted model further used to calculate the modified \emph{AIC} or \emph{BIC}. See references for more details. Note that, for convenience in later usage, returned value is multiplied by -1 inside the function so that large value denotes good quality partition.
}

\value{A numeric value representing the quality of partition under consideration in terms of follow-up.}

\references{
Liang,H. and Zou,G.H. (2008). "Improved AIC selection strategy for survival analysis", \emph{Comput Stat Anal}., 52, 2538-2548.

Volinsky,T.C. and Raftery,A.E. (2000). "Baysian information criteria for censored survival models", \emph{Biometrics}, 56, 256-262.

Obulkasim,A. et al., (2013). "Semi-supervised adaptive-height snipping of the Hierarchical Clustering tree", submitted.
}

\author{Askar Obulkasim}


\seealso{
 \code{\link{measure}}
}

\examples{
data(BullingerLeukemia)
em <- BullingerLeukemia[["em"]]
st <- BullingerLeukemia[["surv.time"]]
event <- BullingerLeukemia[["status"]]
cl <- HCsnipper(em[, 1:60], min = 5)
result <- apply(cl[["partitions"]], 1, function(x) surv_measure(x, st[1:60], event[1:60]))

}


