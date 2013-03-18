\name{perm_test}
\alias{perm_test}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
A function to select the optimal partition (clustering) from large number of candidates and calculate the \emph{p-value} for it. 
}

\description{For a given set of partitions (each partition is composed of non-overlapping clusters), this function uses two types of data to evaluate each partition and select the optimal one which has the highest rank in terms of both data type (presumed that \emph{score1} and \emph{score2} were from two different data source). Permutation approach used to calculate the corrected \emph{p-value} of the selected partition.}

\usage{
perm_test(partitions, surv.time, status, score1 = NULL, score2, method = "BIC", nperm = 1000)
}


\arguments{
  \item{partitions}{A matrix in which rows represent partitions and columns represent samples}
  
  \item{surv.time}{A numeric vector contains follow-up time of patients in the partition}
  
  \item{status}{A binary vector contains survival status of patients in the partition, 0 = alive, 1 = dead}
  
  \item{score1}{A numeric vector contains the quality score for each partition. Scores are assumed to be calculated using the follow-up data. Note, prepare this vector in a way that high value corresponds to good quality partition.}
  
  \item{score2}{A numeric vector contains the quality score for each partition calculated by using any data type except for follow-up. The same as \emph{score1} this vector must be prepared in a way that high value corresponds to good quality partition.}
  
  \item{method}{Type of partition evaluation measure to use. Must be the same as the type of measure used in calculating the \emph{score1}. Default is \emph{BIC}}

  \item{nperm}{The number of permutations.}

}

\details{
When studying association of cluster membership with follow-up data, we cannot use the standard testing procedures. Because \emph{score1} is already used the follow-up data. Thus, we would use the follow-up data twice and the resulting p-value is likely to be too small. We avoid this bias by also applying the semi-supervised partition selection under the null-hypothesis. This null-hypothesis is simply the absence of association between the data type used to generate the 'score2' and the follow-up. Our partition selection in combination with a suitable test statistic is designed to detect associations that can be represented by groups of samples. We adapt the p-value computation as follows:

\enumerate{
\item Use a suitable test statistic (e.g. log-rank for time-to-event data and chi-square for nominal data) to compute the conditional p-value given the cluster labels in the selected partition: \emph{p_obs}.

\item For i = 1...nperm:
    \enumerate{
       \item Randomly permute follow-up data among the samples.
       \item Apply exactly the same type of evaluation measure to evaluate all partitions, e.g. generate new \emph{score1} , but \emph{score2} is fixed. Selected the best
            partition as before.
       \item Conditional on the resulting partition, compute \emph{p-value} \emph{p_i}.
    }    
\item Finally, the \emph{p-value} of interest is equal the number of time \emph{p_obs} smaller than the \emph{p-value}s from the permutation divided by the number of permutations ran.
}

Here, \emph{p} satisfies a crucial property of \emph{p-value}s: it is uniformly distributed when the null-hypothesis is true, because then \emph{p_obs} and \emph{p_i} are exchangeable random variables. The exchangeability is a result from the null-hypothesis and the use of exactly the same procedures to compute \emph{p_obs} and \emph{p_i}. 


}


\value{
A list object contains following objects:
  \item{obs.p}{Observed \emph{p-value}}
 \item{perm.p }{A vector of \emph{p-value}s from permutations.}
 \item{best}{Selected optimal partition}
}

\references{
Obulkasim,A. et al., (2013). "Semi-supervised adaptive-height snipping of the Hierarchical Clustering tree, with application to genomic data", submitted.
}

\author{
Askar Obulkasim
}


\seealso{
\code{\link{TwoHC_perm}}
}

\examples{
data(BullingerLeukemia)
em <- BullingerLeukemia[["em"]]
st <- BullingerLeukemia[["surv.time"]]
event <- BullingerLeukemia[["status"]]
cl <- HCsnipper(em[, 1:60], min = 5)
m <- apply(cl[["partitions"]], 1, function(x) measure(x = x, dis = 1-cor(em[, 1:60]))) 
s <- apply(cl[["partitions"]], 1, function(x) surv_measure(x, st[1:60], event[1:60]))
result <- perm_test(cl[["partitions"]], st[1:60], event[1:60], score1 = s, score2 = m, nperm = 10) 

### Visualize cluster differences in terms of Entropy.
H <- EnvioPlot(dat = em[, 1:60], x= result[["best"]])
} 




