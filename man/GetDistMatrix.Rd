\name{GetDistMatrix}
\alias{GetDistMatrix}
\title{Create a distance matrix}

\usage{
GetDistMatrix(Data, dist_method = "pearson")
}

\description{
Creates a distance matrix from high-dimensional expression data by the chosen distance metric
}

\arguments{
  \item{Data}{
    A matrix or a dataframe of gene expression data
  }
  \item{dist_method}{
    String denoting distance metric ("euclidean","manhattan","pearson" or "spearman"). Default is "pearson".
  }
}

\examples{
  Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
  dist = GetDistMatrix(Data, dist_method = "pearson")
}

