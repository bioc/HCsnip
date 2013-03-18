\name{TcgaGBM}
\alias{TcgaGBM}
\docType{data}

\title{Glioblastoma multiforme gene expression data}

\description{
The subset of latest version of TCGA glioblastoma (GBM) level 3 gene expression data with partial clinical info. Contains expression profiles of 120 samples measured on 3000 genes. Clinical data includes follow-up and type of drugs patients have been administered.}

\usage{data(TcgaGBM)}

\format{An object of list class}


\source{The Cancer Genome Atlas Network (2008), "Comprehensive genomic characterisation defines human glioblastoma genes and core pathways", \emph{Nature}, 490, 61-70.}

\note{The complete dataset is available at the TCGA data portal \url{https://tcga-data.nci.nih.gov/tcga/}.}

\examples{
data(TcgaGBM)
names(TcgaGBM)
}

