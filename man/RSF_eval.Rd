\name{RSF_eval}
\alias{RSF_eval}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Function to calculate error rate using the Random Survival Forest }

\description{This function constructs survival forest using training sample's follow-up as response and cluster labels as covariate. Constructed forest is used to calculate cumulative hazard function (CHF) for each new sample based its cluster label. CHFs are compared with new samples' actual survival time to calculate the error rate. Error rate ranges from 0 to 1, with 0 representing perfect.
}

\usage{
RSF_eval(partition, surv.time, status, te.partition, te.surv.time, te.status, ...)
}


\arguments{
  \item{partition}{Partition (clustering) corresponds to the training set samples. Cluster labels in the partition will be used as covariate to construct SF.}
  
  \item{surv.time}{A numeric vector contains the follow-up information of patients in \emph{partition} , will be used as response to construct SF.}
  
  \item{status}{A binary vector contains survival status of patients in the partition, normally 0=alive, 1=dead, will be used as response to construct SF.}
  
  \item{te.partition}{Partition to be evaluated}
  
  \item{te.surv.time}{A numeric vector contains the follow-up information of patient's in \emph{te.partition}}
  
  \item{te.status}{A binary vector contains survival status of patients in \emph{te.partition}, normally 0=alive, 1=dead.}
  
  \item{\dots}{Arguments for \code{\link{rsf}} from the \pkg{randomSurvivalForest} package}  
}


\value{A vector of error rates. Length is equal the number of trees constructed. Default is 1000.}

\references{
Ishwaran,H. et al., (2008). "Random survival forest", \emph{Ann. App. Statist.}, 2, 841-860.

Obulkasim,A. et al., (2013). "Semi-supervised adaptive-height snipping of the Hierarchical Clustering tree, with application to genomic data", submitted.
}

\author{
Askar Obulkasim
}


\seealso{
\code{\link{rsf}}
}

\examples{
data(BullingerLeukemia)
em <- BullingerLeukemia[["em"]]
st <-  BullingerLeukemia[["surv.time"]]
event <- BullingerLeukemia[["status"]]
cl <- HCsnipper(em[, 1:60], min = 5)
pred <- cluster_pred(X = em[, 1:100], partition = cl$partitions[1, ], surv.time = st[1:60],
                     status = event[1:60], te.index = 61:100)
Err <- RSF_eval(cl$partitions[1, ], st[1:60], event[1:60], pred, st[61:100], event[61:100])

}




