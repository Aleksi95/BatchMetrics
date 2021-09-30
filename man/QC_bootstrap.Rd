\name{QC_bootstrap}
\alias{QC_bootstrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Quality control bootstrap analysis
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
  Compute quality control metrics for batch correction on bootstrap-resampled datasets
}
\usage{
QC_bootstrap(data_list,
            biol.groups,
            batches,
            method_names = NULL,
            metrics = c("F-score",
                        "Davies-Bouldin",
                        "kBET", "kNN",
                        "KL-divergence",
                        "Silhouette",
                        "kNN",
                        "mindist"),
            Fscore_method = "scaled",
            iters = 50,
            savefile = FALSE,
            filename = "evaluations.csv",
            plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_list}{
%%     ~~Describe \code{x} here~~
  A list of datasets of equal dimensions, for example batch-corrected datasets or different levels of batch correction performed on a single dataset (such as produced by the dilute_series function).
}
  \item{biol.groups}{
    A vector containing the cluster labels of the biological variables of interest, the length of columns of the datasets.
  }

  \item{batches}{
    A vector containing the cluster labels of the batch effect variables of interest, the length of columns of the datasets.
  }

  \item{method_names}{
    A vector containing the names of batch correction methods used on the different datasets. Default is method 1, method 2 etc.
  }

  \item{metrics}{
    A vector of strings containing the names of the quality control metrics to run on the datasets. Must be a combination of one or more of the        metrics listed above.
  }

  \item{Fscore_method}{
    (String) A method for computing the f-score metric (see \code{\link{getFscores}}). Must be one of "stratified", "scaled", or "group mean". Default is "scaled".
  }

  \item{iters}{
    A number of iterations for the bootstrap resampling. Default is 50. NOTE: higher number of iterations may be computationally inefficient.
  }

  \item{savefile}{
    Logical. Denotes whether the results should be saved in a csv-format, default is FALSE.
  }

  \item{filename}{
    The name of the csv-file to save.
  }

  \item{plot}{
    Logical. Denotes whether the function automatically produces a boxplot of the bootstrap results. Default is FALSE.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{means }{A dataframe containing the mean for each metric over the bootstrap runs, for biological and batch variables and their ratios}
  \item{results }{The full results of the metrics run on each resampled dataset}
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
#Create a list of datasets to calculate the bootstrapped metrics on

datalist = lapply(1:5, function(d) matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000))
biol.groups = sample(1:4, size = 100, replace = TRUE)
batches = sample(1:2, size = 100, replace = TRUE)

#calculate the bootstrapped metrics
qcmetrics = QC_bootstrap(datalist, biol.groups, batches, plot = FALSE)
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
