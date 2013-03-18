\name{cluster_pred}
\alias{cluster_pred}

\title{Semi-supervised clustering}

\description{For a given partition, this function assigns new samples to one of the clusters in the partition. Partition (clustering) is composed of non-overlapping clusters}

\usage{
cluster_pred(X, partition, surv.time = NULL, status = NULL, te.index, minclus = 4, 
             te.surv.time = NULL, te.status = NULL, method = "conc", maxmiss = 30, 
              plot.it = FALSE, ...)
}

\arguments{
  \item{X}{An object of class \code{\link{ExpressionSet}} or data matrix in which columns are assumed to represent the samples, and rows represents the sample's features. Missing values are allowed. \emph{X} must include the data set from which \emph{partition} is obtained and the data set of test samples.}
  
  \item{partition}{A numeric vector contains the non-overlapping cluster labels.}
  
  \item{surv.time}{A numeric vector contains follow-up time of patients in the partition.}
  
  \item{status}{A binary vector contains survival status of patients in the partition, 0 = alive, 1 = dead.}
  
  \item{te.index}{A numeric vector contains the indices of columns in \emph{X} corresponds to the test samples.}
  
  \item{minclus}{The minimum number samples allowed to form a cluster for the test set. This is to avoid returning tiny clusters and reduce the effect of outliers.}
  
  \item{te.surv.time}{An optional vector contains follow-up time of patients in the test set. If supplied with \emph{te.status}, the logrank test \emph{p-value} is calculated for the test set.}
  
  \item{te.status}{An optional vector contains survival status of patients in the test set. If supplied with \emph{te.surv.time}, the logrank test \emph{p-value} is calculated for the test set.}
  
  \item{method}{Type of methods to use in assigning test samples to one of the clusters in the partition. Must be either Ward distance \emph{ward} or Harrel's concordance index \emph{conc} (default).}

  \item{maxmiss}{Maximum percentage of missing values per row in \emph{X}}
  
  \item{plot.it}{If TRUE and follow-up data of the test samples are given, Kaplan Meier curve(s) will be generated for each cluster in the test set.}
  
  \item{\dots}{Arguments for \code{\link{impute.knn}} from the \pkg{impute} package.}
}  
  
\details{
User has two options to assign test set to one of the clusters in the partition. One option is to use the Ward distance. Specifically, an average distance is calculated between a test sample and samples in each cluster in the partition, separately. The test sample is assigned to a cluster for which average distance is the smallest. Follow-up data is not required for this option. 

Second option is to use the Harrel's concordance index (Harrel et al., 1982). For this option both main and follow-up data corresponds to the given partition are required. Main data is used to find the pseudo nearest neighbours (PNN) of a test sample (Obulkasim et a., 2011), and follow-up data is used to check how much PNN's follow-up info is concordant with follow-up info of samples in each cluster. The test sample is assigned to a cluster for which average concordance is the highest. 


Before selecting either one of the options, we recommend user to check the correlation between main data and follow-up info (e.g. using global test). If correlation is relatively large, we recommend to use 'conc' option, and vice versa.

}

\value{If \emph{plot.it} is FALSE, function returns a vector of cluster labels of the test set. If TRUE and follow-up data of the test set are given, function returns a list object contains following components:
   \item{St}{a data frame with following five columns:
         \tabular{rl}{
           [,1] \tab the unique survival times found in the test set's follow-up.\cr
           [,2] \tab the cluster labels of the test set.\cr
           [,3] \tab the number of patients at risk at each unique time point.\cr
           [,4] \tab the survival probability at each unique time point.\cr
           [,5] \tab the variance of the survival probability at each unique time point.\cr
        }
   }        
   %in which the first column are the unique survival times found in the test set's follow-up; second column is the group labels of the test set; third column is the number of patients at risk at each unique time point; forth column is the survival probability at each unique time point; fifth column is the variance of the survival probability at each unique time point. }
   \item{value}{logrank test \emph{p-value} for the test set}
   
}

\references{
Obulkasim,A. et al., (2013). "Semi-supervised adaptive-height snipping of the Hierarchical Clustering tree, with application to genomic data", submitted.

Harrel,E.F. et al., (1982). "Evaluating the yield of medical tests", \emph{JAMA}, 247, 2543-2546.

Obulkasim,A. et al., (2011). "Stepwise classification of cancer samples using clinical and molecular data", \emph{BMC Bioinformatics}, 12, 422.

Troyanskaya,O. et al., (2001). "Missing value estimation methods for DNA microarrays". \emph{Bioinformatics}, 17, 520-525.
}

\author{Askar Obulkasim}


\seealso{
\code{\link{TwoHC_assign}}
}

\examples{
data(BullingerLeukemia)
em <- BullingerLeukemia[["em"]]
st <-  BullingerLeukemia[["surv.time"]]
event <- BullingerLeukemia[["status"]]
cl <- HCsnipper(em[, 1:60], min = 5)
result <- cluster_pred(X = em[, 1:100], partition = cl$partitions[1, ], surv.time = st[1:60], 
                       status = event[1:60], te.index = 61:100)
names(result)

}









