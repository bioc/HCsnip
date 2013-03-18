\name{EnvioPlot}
\alias{EnvioPlot}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Visualize cluster's molecular entropy by violin plot}

\description{
This function first calculates the entopy of each sample in the given partition, and make a violin plot for each cluster. If clusters are different in term of their molecular profiles, then one may expect desnity differences in violin plots.
}

\usage{
EnvioPlot(dat, method = "knn", x, horizontal = FALSE, 
               col = NULL, names = NULL, ...)
}

\arguments{
  \item{dat}{An object of class \code{\link{ExpressionSet}} or data matrix from which partition is obtained. Columns are assumed to represent samples, and rows represent sample's features. Missing values are allowed.}
  \item{method}{Type of method to calculate sample's molecular entropy. Either \emph{knn} or \emph{normal}.}
  \item{x}{A partition for which violin plot to be made.}
  \item{horizontal}{Should boxes are organized horizontally? default is FALSE.}
  \item{col}{A vector of colors for each cluster. Should be equal the number of clusters in \emph{x}.}
  \item{names}{A vector of labels for each cluster.}
  \item{\dots}{Arguments for \code{\link{hdEntropy}} to calculate the molecular entropy from the \pkg{sigaR} package.}  
}

%\details{The entropy estimate is returned as a numeric. 
%}

\value{The entropy estimate is returned as a numeric.}

\references{
van Wieringen,N.W. and van der Vaart,W.A. (2010) "Statistical analysis of the cancer cell's molecular entropy using high-throughput data", \emph{Bioinformatics}, 27, 556-563.

Adler,D. (2005). vioplot: A violin plot is a combination of a box plot and a kernel density plot, R package, \url{http://cran.r-project.org/web/packages/vioplot/index.html}.
}

\author{Askar Obulkasim}


\seealso{
\code{\link{hdEntropy}}
}

\examples{
data(BullingerLeukemia)
em <- BullingerLeukemia[["em"]]
cl <- cutree(hclust(as.dist(1 - cor(em[, 1:60])), method = "ward"), k = 2)
result <- EnvioPlot(dat = em[, 1:60], x = cl)
}
