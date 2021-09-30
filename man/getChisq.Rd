\name{getChisq}
\alias{getChisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Chi-square test for clustering accuracy.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A metric for cluster separation that uses a hierarchical clustering algorithm and compares clustering results with chi-square contingency table test.
}
\usage{
getChisq(dists, sample_types)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{dists}{
%%     ~~Describe \code{x} here~~
  A distance matrix or a list of distance matrices, such as one produced by GetDistMatrix
}
  \item{sample_types}{
  A sample type vector containing the cluster labels for the sample types of interest (biological groups or batches)
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{statistics}{a value or a vector of value of the chi-square test statistics}
\item{pvals}{a value or a vector of value of the p-values of the chi-square test}
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

#distance matrices
dists = lapply(datalist, function(d) GetDistMatrix(d, dist_method = "pearson"))
sample_types = sample(1:2, 100, replace = TRUE)

getChisq(dists, sample_types)

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
