\name{DaviesBouldinScores}
\alias{DaviesBouldinScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Davies-Bouldin indices for many datasets
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A wrapper of \code{\link{DaviesBouldinIndex}} for one or multiple datasets
}
\usage{
DaviesBouldinScores(dists, sample_types)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dists}{
%%     ~~Describe \code{x} here~~
  A distance matrix or a list of distance matrices
}

  \item{sample_types}{
    A sample type vector containing the cluster labels for the sample types of interest (biological groups or batches).
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
%% ...
A vector of Davies-Bouldin indices for the datasets
}
\references{
%% ~put references to the literature/web site here ~
D. L. Davies and D. W. Bouldin, "A Cluster Separation Measure," in IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. PAMI-1, no. 2, pp. 224-227, April 1979, doi: 10.1109/TPAMI.1979.4766909.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~
\code{\link{DaviesBouldinIndex}}
}
\examples{
##Create a list of datasets to calculate the metric on

datalist = lapply(1:5, function(d) matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000))
#distance matrices
dists = lapply(datalist, function(d) GetDistMatrix(d, dist_method = "pearson"))
sample_types = sample(1:2, 100, replace = TRUE)

DaviesBouldinScores(dists, sample_types)
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
