\name{getkBET}
\alias{getkBET}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
kBET rejection rate
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes kBET (BÃ¼ttner et al 2019) rejection rate and average p-value for one or multiple datasets.
}
\usage{
getkBET(dists, sample_types)
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
%%  \item{comp1 }{Description of 'comp1'}
\item{avg.pvals}{The kBET average p-value (see \code{\link{kBET}})}
%%  \item{comp2 }{Description of 'comp2'}
\item{rej.rates}{The rejection rates of kBET (see \code{\link{kBET}})}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
BÃ¼ttner, M., Miao, Z., Wolf, F.A. et al. A test metric for assessing single-cell RNA-seq batch correction. Nat Methods 16, 43-49 (2019). https://doi.org/10.1038/s41592-018-0254-1
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

require(FNN)
require(kBET)

#Create a list of datasets to calculate the metric on

datalist = lapply(1:5, function(d) matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000))

#distance matrices
dists = lapply(datalist, function(d) GetDistMatrix(d, dist_method = "pearson"))
sample_types = sample(1:2, 100, replace = TRUE)

getkBET(dists, sample_types)
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
