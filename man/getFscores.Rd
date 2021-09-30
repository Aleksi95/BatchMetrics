\name{getFscores}
\alias{getFscores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
F-scores for multiple datasets
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A wrapper of \code{\link{getFscore}} for multiple datasets.
}
\usage{
getFscores(dists, sample_types, method = "scaled", quantile = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dists}{
%%     ~~Describe \code{x} here~~
  A list of distance matrices
}

  \item{sample_types}{
    A sample type vector or a list of vectors containing the cluster labels for the sample types of interest (biological groups or batches).
  }

  \item{method}{
    The method of computing the f-score (see \code{\link{getFscore}}). one of ("scaled", "stratified" or "group mean")
  }

  \item{quantile}{
    The quantile or a vector of quantiles to scale the distances by
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
A vector containing the f-scores for multiple datasets
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
\code{\link{getFscore}}
}
\examples{
##Create a list of datasets to calculate the metric on

datalist = lapply(1:5, function(d) matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000))
dists = lapply(datalist, function(d) GetDistMatrix(d, dist_method = "pearson"))
sample_types = sample(1:2, 100, replace = TRUE)

getFscores(dists, sample_types)

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
