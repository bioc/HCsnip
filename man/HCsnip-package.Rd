\name{HCsnip-package}
\alias{HCsnip-package}
\docType{package}
\title{Semi-supervised adaptive-height snipping of the Hierarchical Clustering tree}

\description{This package contain functions for extracting meaningful clusters from a HC tree. Rather than cutting the HC tree at a fixed highest (as existing methods do), it snips the tree at variable heights to extract hidden clusters. Cluster extraction process uses both the data type from which HC tree is derived and the available patients follow-up information. Functions for testing the significance of extracted clusters and cluster visualization using sample's molecular entropy are also given. If two HC trees are presented which maybe corresponding to the two treatment groups, this package also includes functions for optimally assigning new patients to one of the two HC trees and calculate the expected gain in terms of follow-up.}    


\details{
\tabular{ll}{
Package: \tab HCsnip\cr
Type: \tab Package\cr
License: \tab GPL (>= 2)\cr
}

}

\author{ Askar Obulkasim
Maintainer: Askar Obulkasim <askar.wubulikasimu@vumc.nl>

}

\references{
Obulkasim,A. et al., (2013). "Semi-supervised adaptive-height snipping of the Hierarchical Clustering tree", submitted.
}
