\name{TwoHC_assign}
\alias{TwoHC_assign}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Function to assign new samples to one of the two given hierarchical clustering (HC) trees in a semi-supervised way}

\description{For given molecular data sets from two non-overlapping groups of patients, this functions constructs two independent HC trees and assigns new samples to one of them in semi-supervised way. See details.}

\usage{
TwoHC_assign(X, index1, index2, new.X, dis.method = "cor", link.method = "ward", 
             minclus = 4, maxmiss = 30, surv.time, status, method1 = "BIC", 
              method2 = "g2")
}

\arguments{
  \item{X}{An object of class \code{\link{ExpressionSet}} or data matrix from which two HC tress to be derived. Columns are assumed to represent the samples, and rows represent the sample's features. Missing values are allowed.}
  
  \item{index1}{Column indices of patients in \emph{X} correspond to the first group.}
  
  \item{index2}{Column indices of patients in \emph{X} correspond to the second group.}
  
  \item{new.X}{An object of class \code{\link{ExpressionSet}} or data matrix corresponds to new samples. Columns are assumed to represent the samples, and rows represents the sample's features. Missing values are allowed.}
  
  \item{dis.method}{The distance measure to be used. This must be one of method acceptable for \code{\link{dist}} function or the Pearson correlation (default).}
  
  \item{link.method}{The agglomeration method to be used. This should be one of "ward" (default), "single", "complete", "average", "mcquitty", "median" or "centroid".}
  
  \item{minclus}{The minimum number of samples allowed to form a cluster. This parameter inversely proportional to the number of partition returned from a HC tree. e.g. a large value returns small number of partitions, and vice versa.}
  
  \item{maxmiss}{Maximum percentage of missing values per row in \emph{X}.}
  
  \item{surv.time}{A numeric vector contains follow-up information of patient's in \emph{X} }
  
  \item{status}{A binary vector contains survival status of patients in \emph{X}, normally 0=alive, 1=dead.}
  
  \item{method1}{Type of partition evaluation measures to use for assessing the relationship between follow-up and a partition. Default is "BIC".}
  
  \item{method2}{Type of Partition evaluation measure to use for assessing the relationship between data matrix \emph{X} and a partition. Default is Goodman and Kruskal index "g2".}
  
}

\details{
Say molecular profiles of two groups patients (without overlap) treated with two different drugs or the same drugs in different combinations are available. Besides that, their follow-up information are also given. When a new patient comes in, question will be to which group this patient should be assigned so that he/she will benefit most by the type of treatment this group received. 

This function is designed for this problem. it works as follows: first, two independent HC trees will be derived from given data; second, partitions are extracted and the optimal partition is selected from each HC tree, separately; third, new patient's molecular profile (presumed that new patient's follow-up and group label are missing) is compared with each cluster in each optimal partition to calculate average similarity and identify two most similar clusters (competing clusters); finally, new sample is assigned to one of the two competing clusters which has better survival.

}

\value{A list object contains following components:
  \item{hc1}{HC tree derived from the data corresponds to the first group.}
  \item{hc2}{HC tree derived from the data corresponds to the second group.}
  \item{partitions.hc1}{A matrix includes the partition extracted from \emph{hc1}. Rows represent partitions and columns represent samples.}
  \item{partitions.hc2}{A matrix includes the partition extracted from \emph{hc2}. Rows represent partitions and columns represent samples.}
  \item{best.hc1}{Optimal partition selected from \emph{partitions.hc1}} 
  \item{best.hc2}{Optimal partition selected from \emph{partitions.hc2}} 
  \item{score.hc1}{A matrix with two columns. The first column contains the quality scores of \emph{partitions.hc1} calculated using the follow-up data. The second column contains the quality scores of \emph{partition.hc1} calculated by using \emph{X}.}
  \item{score.hc2}{The same as \emph{score.hc1}, but for \emph{partitions.hc2}.}
  \item{Assign}{A matrix with three columns. The first column contains the indices of HC trees to which a test sample was assigned. The second column contains the indices of clusters in \emph{best.hc1} to which a test sample was most similar. The third column contains the indices of clusters in \emph{best.hc2} to which a test sample was most similar.}
  \item{surv.time}{The same as input}
  \item{status}{The same as input}
  \item{index1}{The same as input}
  \item{index2}{The same as input}
  \item{new.X}{The same as input}
  \item{X}{The same as input}
  \item{method1}{The same as input}
  \item{method2}{The same as input}
  \item{minclus}{The same as input}
   
}

\references{
Harrel,E.F. et al., (1982). "Evaluating the yield of medical tests", \emph{JAMA}, 247, 2543-2546.

Obulkasim,A. et al., (2011). "Stepwise classification of cancer samples using clinical and molecular data", \emph{BMC Bioinformatics}, 12, 422.

Troyanskaya,O. et al., (2001). "Missing value estimation methods for DNA microarrays". \emph{Bioinformatics}, 17, 520-525.

Obulkasim,A. et al., (2013). "Semi-supervised adaptive-height snipping of the Hierarchical Clustering tree, with application to genomic data", submitted.
}

\author{Askar Obulkasim}


\seealso{
See also \code{\link{TwoHC_perm}}, \code{\link{cluster_pred}}
}

\examples{
data(TcgaGBM)
id1 <- which(TcgaGBM[["drugs"]] == "Avastin")
id2 <- which(TcgaGBM[["drugs"]] == "Temodar") 
em <- TcgaGBM[["em"]]
st <-  TcgaGBM[["surv.time"]]
event <- TcgaGBM[["status"]]
result <- TwoHC_assign(X = em[ ,c(id1[1:50], id2[1:50])], index1 = 1:50, index2 = 51:100, 
                      new.X = em[, c(id1[51:60], id2[51:60])], minclus = 4,
                     surv.time = st[c(id1[1:50], id2[1:50])], 
                     status = event[c(id1[1:50], id2[1:50])])  
}
