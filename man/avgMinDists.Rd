\name{avgMinDist}
\alias{avgMinDist}
\alias{getMinDists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Average minimum cluster distance
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A metric for cluster separation/batch effect evaluation that computes the average minimum distance between batches
}
\usage{
avgMinDist(dists, batch)
getMinDists(dists, batch)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{dists}{
%%     ~~Describe \code{x} here~~
  A distance matrix or a list of distance matrices (in case of getMinDists)
}

  \item{batch}{
    A vector of labels (factor/numeric/character) describing the batch effects.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
getMinDists is a wrapper for multiple datasets
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A value or a vector of values of the average minimum distance between the batches.
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
#A list of datasets
datalist = lapply(1:5, function(d) matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000))
#distance matrices
dists = lapply(datalist, function(d) GetDistMatrix(d, dist_method = "pearson"))
#batch labels
batch = sample(1:2, 100, replace = TRUE)

#for a single dataset
avgMinDist(dists[[1]], batch)

#for multiple datasets
getMinDists(dists, batch)
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
