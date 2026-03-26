\name{evalBatchEffect}
\alias{evalBatchEffect}
\alias{runMetrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluate batch effects with a chosen metric
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{evalBatchEffect} computes a chosen cluster separation metric on a dataset for both biological groups and batch labels, and returns the metric values and their ratio. \code{runMetrics} is a wrapper that runs multiple metrics at once.
}
\usage{
evalBatchEffect(data, sample_types, batch, metric = "Davies-Bouldin",
                y = NULL, distMatrix = NULL, zeroRows = FALSE)

runMetrics(data, sample_types, batch, y = NULL, distMatrix = NULL,
           metrics, zeroRows = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
A gene expression data matrix with rows as features (genes) and columns as samples.
}

  \item{sample_types}{
A vector of labels for the biological sample groups.
}

  \item{batch}{
A vector of batch labels for each sample.
}

  \item{metric}{
String specifying the metric to compute. One of \code{"F-score"}, \code{"Davies-Bouldin"}, \code{"kBET"}, \code{"Silhouette"}, \code{"kNN"}, \code{"mindist"}, \code{"kldist"} or \code{"gPCA"}. Default is \code{"Davies-Bouldin"}.
}

  \item{y}{
Optional response variable. Currently used only for the \code{"pvca"} metric.
}

  \item{distMatrix}{
An optional pre-computed distance matrix. If NULL, it will be computed from \code{data} using \code{\link{GetDistMatrix}}.
}

  \item{zeroRows}{
Logical. Whether to include rows with zero variance. Default is FALSE.
}

  \item{metrics}{
A character vector of metric names to compute (used by \code{runMetrics}).
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{evalBatchEffect} computes the chosen metric for both the biological groups (\code{sample_types}) and batch labels, and returns a named numeric vector with three elements: the metric value for biological groups (\code{bio}), the metric value for batches (\code{batch}), and their ratio (\code{ratio}). A higher ratio indicates better separation of biological groups relative to batches.

\code{runMetrics} calls \code{evalBatchEffect} for each metric in the \code{metrics} vector and returns the results as a named list.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
For \code{evalBatchEffect}: A named numeric vector with elements \code{bio}, \code{batch} and \code{ratio}.

For \code{runMetrics}: A named list where each element is the result of \code{evalBatchEffect} for the corresponding metric.
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
\code{\link{QC_bootstrap}}, \code{\link{getFscore}}, \code{\link{DaviesBouldinIndex}}, \code{\link{getSilhouette}}
}
\examples{
#Generate a gene expression data matrix (genes x samples)
Data = matrix(rnorm(500 * 30, mean = 0, sd = 1), nrow = 500, ncol = 30)

#Define sample type and batch labels
sample_types = sample(c("A", "B"), 30, replace = TRUE)
batch = sample(c("batch1", "batch2"), 30, replace = TRUE)

#Evaluate one metric
result = evalBatchEffect(Data, sample_types, batch, metric = "F-score")
result

#Evaluate multiple metrics using runMetrics
results = runMetrics(Data, sample_types, batch,
                     metrics = c("F-score", "Davies-Bouldin"))
results
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
