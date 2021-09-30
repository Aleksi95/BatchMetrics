\name{DaviesBouldinIndex}
\alias{DaviesBouldinIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Davies-Bouldin index of cluster separation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Returns the Davies-Bouldin index of cluster separation. Higher values indicate better cluster separation
}
\usage{
DaviesBouldinIndex(dist, sample_types)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{
%%     ~~Describe \code{x} here~~
  A distance matrix, such as one produced by GetDistMatrix
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
%%  \item{comp2 }{Description of 'comp2'}
  \item{DB_score}{The Davies Bouldin score defined as 1 over the mean ratio of between and within cluster distances}

  \item{DB_score2}{An alternative way of computing the Davies-Bouldin score}
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Generate dataset and sample types

Data = matrix(rnorm(100*10,mean=0,sd=1), ncol=10, nrow=100)
dist = GetDistMatrix(Data, dist_method = "pearson")
sample_types = sample(1:2, 10, replace = TRUE)

#Calculate the metric
DaviesBouldinIndex(dist, sample_types)
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
