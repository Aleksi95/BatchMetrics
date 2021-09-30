\name{getSilhouettes}
\alias{getSilhouettes}
\alias{getSilhouette}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Silhouette score for cluster separation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates the silhouette cluster separation score given a distance matrix (or matrices) of a dataset or multiple datasets
}
\usage{
getSilhouettes(dist, sample_types)
getSilhouette(dist, sample_types)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{
%%     ~~Describe \code{x} here~~
A distance matrix or a list of distance matrices (for getSilhouettes).
}

\item{sample_types}{
A vector of labels specifying the cluster labels (sample types or batches).
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
This metric of cluster separation is similar to \code{\link{DaviesBouldinIndex}}. The silhouette value of a single data point is defined as

\deqn{
s(i) = \frac{a(i) - b(i)}{max\{a(i), b(i)\}}
}

for clusters with more than one datapoint an \eqn{s(i)} for clusters with one data point (a(i) is the mean distances between data points within a cluster and b(i) is the minimum mean distance to points in other clusters).
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
For getSilhouette:

\item{MeanSilh}{Mean silhouette values}
\item{SilhouetteVector}{The vetor of all silhouette values}
%% ...

For getSilhouettes:
The vector of mean silhouette values (MeanSilh) for all given datasets.
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
\code{\link{DaviesBouldinIndex}}
}
\examples{
#For a single dataset

Data = matrix(rnorm(100*10,mean=0,sd=1), ncol=10, nrow=100)
dist = GetDistMatrix(Data, dist_method = "pearson")
sample_types = sample(1:2, 10, replace = TRUE)

getSilhouette(dist, sample_types)

###for multiple datasets

datalist = lapply(1:5, function(d) matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000))
dists = lapply(datalist, function(d) GetDistMatrix(d, dist_method = "pearson"))

getSilhouettes(dist, sample_types)
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
