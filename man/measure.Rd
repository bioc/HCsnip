\name{measure}
\alias{measure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate cluster quality 
}

\description{
Function to evaluate the overall quality of a parition (composed of non-overlapping clusters) by user defined criteria.
}

\usage{
measure(parti, dis, X = NULL, method = "g2", maxmiss = 30, ...)
}

\arguments{

  \item{parti}{Partition to be evaluated. }

  \item{dis}{A square distance matrix or class object of \code{dist} corresponding to \emph{x}.}
  
  \item{X}{data matrix corresponding to the \emph{parti}. Columns are assumed to represent the samples, and rows represent the sample's features. Missing values are allowed. This is an optional argument, but If type is set to 'igp', then matrix must be given.}
  
  \item{method}{Type of evaluation measure to use for assessing the quality of clusters in \emph{x}. Default is Goodman and Kruskal index \emph{g2}.}
  
  \item{maxmiss}{Maximum percentage of missing values per row in \emph{dat}}
  
  \item{\dots}{Arguments for function \code{\link{cluster.stats}} from the \pkg{fpc} package. See details below.}
}

\details{
Numerous cluster quality measuring criteria have been proposed. This package includes only a few well known ones. Except for the 'c.index' and the in group proportion 'igp', rest of the criteria come from the function \code{\link{cluster.stats}} in \pkg{fpc} package. For latter one, please see the returned arguments of the \code{\link{cluster.stats}} function before you decide which criteria to choose. Note that, the value returned by different criteria has different meaning. For example. the larger the Goodman and Kruskal index 'g2' the better, for the index G3 'g3' the smaller the better. Thus, interpret returned value accordingly.
}

\value{
A numeric value representing the quality of partition under consideration.
}

\references{
Hennig,C. (2010). fpc: Flexible procedures for clustering, R package, \url{http://CRAN.R-project.org/package=fpc}.

Kapp,A.V. and Tibshirani,R. (2007) "Are clusters found in one dataset present in another dataset?", \emph{Biostatistics}, 8, 9-31.

Obulkasim,A. et al., (2013). "Semi-supervised adaptive-height snipping of the Hierarchical Clustering tree", submitted.
}


\author{Askar Obulkasim}



\seealso{
\code{\link{surv_measure}}
}

\examples{

data(BullingerLeukemia)
attach(BullingerLeukemia)
cl <- HCsnipper(em[, 1:30], minclus = 5)
cl <- cl$partitions[cl$id, ]
m <- apply(cl, 1, function(x) measure(parti = x, dis = 1 - cor(em[, 1:30]))) 

}



