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
QC_bootstrap(DataNorm,
            biol.groups,
            batches,
            method_names = NULL,
            metrics = c("F-score", "Davies-Bouldin", "kBET", "kNN"),
            levels = seq(0, 1, by = 0.1),
            dist_method = "pearson",
            scaledF = FALSE,
            iters = 50,
            corrMethod = "ComBat",
            dilute_samples = FALSE,
            parallel = TRUE,
            nCores = 16,
            zeroRows = FALSE,
            usePCA = FALSE,
            nPCs = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataNorm}{
  A normalized data matrix (rows are features/genes, columns are samples) to compute bootstrap quality control metrics on.
}
  \item{biol.groups}{
    A vector containing the cluster labels of the biological variables of interest, the length of columns of the dataset.
  }

  \item{batches}{
    A vector containing the cluster labels of the batch effect variables of interest, the length of columns of the dataset.
  }

  \item{method_names}{
    A vector containing the names of batch correction methods. Default is NULL.
  }

  \item{metrics}{
    A vector of strings containing the names of the quality control metrics to run on the dataset. Must be a combination of one or more of "F-score", "Davies-Bouldin", "kBET", "kNN", "mindist", "kldist", "gPCA", "CMS".
  }

  \item{levels}{
    A numeric vector of batch correction levels (proportions between 0 and 1) to evaluate. Default is \code{seq(0, 1, by = 0.1)}.
  }

  \item{dist_method}{
    String denoting the distance method to use for computing the distance matrix. Default is "pearson".
  }

  \item{scaledF}{
    Logical. Whether to use scaled F-score. Default is FALSE.
  }

  \item{iters}{
    A number of iterations for the bootstrap resampling. Default is 50. NOTE: higher number of iterations may be computationally inefficient.
  }

  \item{corrMethod}{
    String denoting the batch correction method to apply. One of "ComBat", "limma", or "harmony". Default is "ComBat".
  }

  \item{dilute_samples}{
    Logical. Whether to use sample dilution series instead of gene bootstrapping. Default is FALSE.
  }

  \item{parallel}{
    Logical. Whether to run bootstrap iterations in parallel. Default is TRUE.
  }

  \item{nCores}{
    Integer. Number of parallel cores to use when \code{parallel = TRUE}. Default is 16.
  }

  \item{zeroRows}{
    Logical. Whether to keep rows with zero standard deviation. Default is FALSE.
  }

  \item{usePCA}{
    Logical. Whether to reduce dimensionality with PCA before computing metrics. Default is FALSE.
  }

  \item{nPCs}{
    Integer. Number of principal components to use when \code{usePCA = TRUE}. Default is 50.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{result}{A list of bootstrap iteration results, each containing bio, batch, and ratio sub-lists for the evaluated metrics at each correction level.}
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
#Create a normalized data matrix to calculate the bootstrapped metrics on

DataNorm = matrix(rnorm(10*100,mean=0,sd=1), ncol=10, nrow=100)
biol.groups = sample(1:2, size = 10, replace = TRUE)
batches = sample(1:2, size = 10, replace = TRUE)

#calculate the bootstrapped metrics (single iteration, not parallel, for illustration)
\dontrun{
qcmetrics = QC_bootstrap(DataNorm, biol.groups, batches, iters = 1, parallel = FALSE)
}
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
