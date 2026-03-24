\name{GetDistMatrix}
\alias{GetDistMatrix}
\title{Create a distance matrix}

\usage{
GetDistMatrix(Data, dist_method = "pearson", fastCor = TRUE)
}

\description{
Creates a distance matrix from high-dimensional expression data by the chosen distance metric
}

\arguments{
  \item{Data}{
    A matrix or a dataframe of gene expression data (rows are features/genes, columns are samples). Must not contain NA values.
  }
  \item{dist_method}{
    String denoting distance metric ("pearson" or "spearman"). Default is "pearson".
  }
  \item{fastCor}{
    Logical. Whether to use the faster WGCNA correlation implementation when available. Default is TRUE.
  }
}

\examples{
  Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
  dist = GetDistMatrix(Data, dist_method = "pearson")
}

