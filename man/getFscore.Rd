\name{getFscore}
\alias{getFscore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
F-score of cluster separation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Computes the f-score for cluster separation by two different methods. Higher f-score indicates better cluster separation
}
\usage{
getFscore(dist, sample_types, method = "scaled", quantile = 0.5)
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

  \item{method}{
    String indicating method of computation ("scaled" or "stratified"). Default is "scaled".
  }

  \item{quantile}{
    The quantile to scale the pairwise distances by. Default is 0.5.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
  The "stratified" method computes the f-score as a product of mean distance ratios of within/between cluster distances for each cluster:

  \deqn{
  F = \prod_{k=1}^K \frac{1/M(\sum_{k \neq n}^M D_{k,n})}{1/M(\sum_{k}^M D_{k,k}}
  }

  where K is the total number of clusters, M is the size of each cluster and D is the distance matrix.

  The "scaled" method scales each pairwise distances by achosen quantile.

  \deqn{
  F = \frac{\sum D_{k, j}/Q(D_k)}{\sum D_{k, k}/Q(D_k)}
  }


}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
The f-score for cluster separation, lower scores indicate better cluster mixing (for example measuring batch effects) and higher score indicates better cluster separation
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
\code{\link{getFscores}}
}

\examples{
#Generate data and sample type vector
Data = matrix(rnorm(100*10,mean=0,sd=1), ncol=10, nrow=100)
dist = GetDistMatrix(Data, dist_method = "pearson")
sample_types = sample(1:2, 10, replace = TRUE)


#Calculate the metric
getFscore(dist, sample_types, method = "scaled", quantile = 0.5)
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
