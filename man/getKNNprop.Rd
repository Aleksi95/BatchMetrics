\name{getKNNprop}
\alias{getKNNprop}
\alias{kNN_proportions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
K-nearest neighbour proportion
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A cluster separation metric that calculates the (mean) proportion of k-nearest neighbours for each sample/datapoint that belong to the same batch or biological cluster of interest.
}
\usage{
getKNNprop(dist, sample_types)
kNN_proportions(dist, sample_types)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{dist}{
%%     ~~Describe \code{x} here~~
  A distance matrix or a list of distance matrices (in case of kNN_proportions)
}

  \item{sample_types}{
    A sample type vector containing the cluster labels for the sample types of interest (biological groups or batches).
  }

}

\details{
%%  ~~ If necessary, more details than the description above ~~
kNN_proportions is a wrapper for multiple datasets
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A value or a vector of values of the mean proportion of datapoints of the same batch/biological group in the k-nearest neighbours. Neighbourhood size is set as the average cluster size of the batches/biological groups.
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##Create a list of datasets to calculate the metric on

datalist = lapply(1:5, function(d) matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000))
dists = lapply(datalist, function(d) GetDistMatrix(d, dist_method = "pearson"))

#sample type vector
sample_types = sample(1:2, 100, replace = TRUE)

#For multiple datasets
kNN_proportions(dists, sample_types)

#For a single dataset
getKNNprop(dists[[1]], sample_types)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
