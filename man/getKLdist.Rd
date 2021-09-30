\name{getKLdist}
\alias{getKLdist}
\alias{getKLdistances}
\alias{getKLD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Kullback-Leibler distance on multiple datasets
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A cluster separation metric that utilises the Kullback-Leibler symmetric distance (FNN package) on one or more datasets
}
\usage{
getKLdist(dists, sample_types, k = 5)
getKLdistances(dists, sample_types, k = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{dists}{
%%     ~~Describe \code{x} here~~
  A distance matrix or a list of distance matrices (in case of getKLdistances)
}

  \item{sample_types}{
    A sample type vector containing the cluster labels for the sample types of interest (biological groups or batches).
  }

  \item{k}{
    The maximum number of nearst neighbours (see KLdist function of the FNN package). The default is 5.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A value or a vector of values of the mean Kullback-Leibler distance between the batches
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
getKLdistances(dists, sample_types, k = 5)
#For a single dataset
getKLdist(dists[[1]], sample_types, k = 5)
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
