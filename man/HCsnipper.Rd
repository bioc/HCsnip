\name{HCsnipper}
\alias{HCsnipper}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{HC tree snipper}

\description{
This function snips given hierarchical clustering (HC) at variable heights to extract all possible partitions. Each partition (clustering) is composed of non-overlapping clusters.
}

\usage{
HCsnipper(X, hc = NULL, dis = NULL, dis.method = "cor", link.method = "ward", 
          minclus = 4, maxmiss = 30, ...)
}

\arguments{
  \item{X}{An object of class \code{\link{ExpressionSet}} or data matrix from which HC tree will be derived. Columns are assumed to represent the samples, and rows represent the sample's features (genes). Missing values are allowed.}
  \item{hc}{HC tree from which partitions to be extracted. Must be an object class of \code{\link{hclust}}. This is an optional argument, but if given \emph{X} and \emph{dis} will be ignored.}
  
  \item{dis}{A square distance matrix or object class of \code{\link{dist}} from which HC tree to be derived. This is an optional argument, if given \emph{X} will be ignored.}

  \item{dis.method}{The distance measure to be used. This must be one of methods acceptable for \code{\link{dist}} function or the Pearson correlation (default).}

  \item{link.method}{The agglomeration method to be used. This should be one of "ward" (default), "single", "complete", "average", "mcquitty", "median" or "centroid".}

  \item{minclus}{The minimum number of samples allowed to form a cluster. This parameter is inversely proportional to the number of partitions returned. e.g. large values returns small number clusters, and vice versa.}

  \item{maxmiss}{Maximum percentage of missing values per row in \emph{X}}
  
  \item{\dots}{Arguments for \code{\link{impute.knn}} from the impute package for missing values imputation in \emph{X}.}
  
}

\details{
For given HC tree, this function snips it at all possible places to extract partitions under the following conditions: 
\itemize{
  \item Singleton is not allowed. 
  \item Snipping places are chosen so that only the samples which are neighbours in the leaf node ordering (see order(hc)) are allowed to form a cluster. 
}  
The last constraint guarantees that sniping does not change the HC tree structure considerably. For example, samples located in far left in the HC tree will not be joined with samples located in far right. The number of partitions return by function depends not only on the \emph{minclus} argument, but also the shape of the HC tree. Large number of partitions can be returned from a balanced HC tree than a skewed one.
}

\value{
This function returns a list object contains following objects:
  \item{partitions}{a matrix in which rows represent partitions and columns represent samples.}
  \item{hc}{HC tree from which partitions are extracted.}
  \item{dat}{data matrix. If \emph{X} has missing values, this will be missing values imputed full data matrix.}
  \item{dis}{the distance matrix used}
  \item{dis.m}{the distance measure used}
  \item{link.m}{the agglomeration method used}

}

\references{
Obulkasim,A. et al., (2013). "Semi-supervised adaptive-height snipping of the Hierarchical Clustering tree, with application to genomic data", submitted.

Troyanskaya,O. et al., (2001). "Missing value estimation methods for DNA microarrays". \emph{Bioinformatics}, 17, 520-525.
}

\author{
Askar Obulkasim
}


\examples{
data(BullingerLeukemia)
em <- BullingerLeukemia[["em"]]
cl <- HCsnipper(em[, 1:60], minclus = 5)
## To check returned partitions size
table(apply(cl[["partitions"]], 1, function(x) length(unique(x))))  
}



