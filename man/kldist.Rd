\name{kldist}
\alias{kldist}
\alias{kldistTwo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kullback-Leibler distance between batches
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{kldist} computes the weighted mean Kullback-Leibler (KL) distance between batch groups across all pairwise batch combinations. \code{kldistTwo} is a helper function that computes the KL distance between two specific batch groups.
}
\usage{
kldist(xba, batch)

kldistTwo(xb1, xb2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xba}{
%%     ~~Describe \code{x} here~~
A numeric matrix with rows as samples and columns as features. The data is scaled before KL distance computation.
}

  \item{batch}{
A factor vector of batch labels for each sample.
}

  \item{xb1}{
A numeric matrix of samples from the first batch (rows as samples, columns as features).
}

  \item{xb2}{
A numeric matrix of samples from the second batch (rows as samples, columns as features).
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{kldist} first standardises the data matrix, then computes pairwise KL distances for all batch pairs using \code{kldistTwo}. The overall result is the weighted mean of these pairwise distances, weighted by the combined sample sizes of each pair, normalised by the total number of samples and the number of batch comparisons.

\code{kldistTwo} computes the KL distance between two groups by comparing the distribution of within-group pairwise distances to the distribution of between-group distances using \code{FNN::KL.dist}. The result is a weighted average of the two directional KL distances.

Note: for \code{kldist}, the \code{batch} argument must be a factor. See \code{\link{getKLdist}} for a version that works with distance matrices.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A single numeric value representing:
\itemize{
  \item For \code{kldist}: the weighted mean KL distance between all batch pairs.
  \item For \code{kldistTwo}: the KL distance between the two batch groups.
}
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
\code{\link{getKLdist}}, \code{\link{getKLdistances}}
}
\examples{
#Generate a data matrix (samples x features, as required by kldist)
Data = matrix(rnorm(20 * 100, mean = 0, sd = 1), nrow = 20, ncol = 100)

#Define batch labels as a factor
batch = factor(sample(c("batch1", "batch2"), 20, replace = TRUE))

#Compute KL distance across all batch pairs
kl = kldist(Data, batch)
kl

#Compute the KL distance between two specific groups
batch_idx1 = which(batch == levels(batch)[1])
batch_idx2 = which(batch == levels(batch)[2])
kl_two = kldistTwo(Data[batch_idx1, ], Data[batch_idx2, ])
kl_two
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
