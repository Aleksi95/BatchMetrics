\name{QC_wrapper}
\alias{QC_wrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quality control wrapper for batch correction evaluation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A comprehensive wrapper function that normalises RNA-seq count data, applies batch correction, and evaluates the quality of the batch correction using multiple metrics. Supports both bootstrapping and resampling-based variance estimation.
}
\usage{
QC_wrapper(CountData, batch, group, y = NULL,
           metrics = c("F-score", "Davies-Bouldin", "kNN", "mindist",
                       "kldist", "gPCA"),
           iters = 100,
           perc = seq(0, 1, by = 0.1),
           corrData = NULL,
           corrMethod = "ComBat",
           deseqData = NULL,
           parallel = TRUE,
           var_measure = c("bootstrap", "resampling",
                           "bootstrapSamples", "resampleSamples"),
           cores = 4,
           design = NULL,
           Ctrl_sample = "normal",
           rnaSeq = "bulk",
           usePCA = TRUE,
           nPCs = 50,
           var.model = "multinom")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CountData}{
%%     ~~Describe \code{x} here~~
A raw RNA-seq count matrix with rows as genes and columns as samples.
}

  \item{batch}{
A vector of batch labels for each sample.
}

  \item{group}{
A vector of biological group labels for each sample.
}

  \item{y}{
An optional response variable. Defaults to \code{group} if NULL.
}

  \item{metrics}{
A character vector of metric names to compute. Available metrics include \code{"F-score"}, \code{"Davies-Bouldin"}, \code{"kNN"}, \code{"mindist"}, \code{"kldist"} and \code{"gPCA"}. Default is all six.
}

  \item{iters}{
The number of bootstrap or resampling iterations. Default is 100.
}

  \item{perc}{
The vector of batch correction percentages for the dilute series. Default is \code{seq(0, 1, by = 0.1)}.
}

  \item{corrData}{
An optional pre-computed batch-corrected data matrix. If NULL, batch correction is performed internally.
}

  \item{corrMethod}{
String specifying the batch correction method to use (\code{"ComBat"}, \code{"limma"} or \code{"harmony"}). Default is \code{"ComBat"}.
}

  \item{deseqData}{
An optional pre-computed DESeq2 object. If NULL and \code{rnaSeq = "bulk"}, DESeq2 normalisation is run internally.
}

  \item{parallel}{
Logical. Whether to run iterations in parallel using \code{mclapply}. Default is TRUE.
}

  \item{var_measure}{
A character vector specifying which variance estimation methods to use. One or more of \code{"bootstrap"}, \code{"resampling"}, \code{"bootstrapSamples"} and \code{"resampleSamples"}. Default is all four.
}

  \item{cores}{
The number of CPU cores to use for parallel computation. Default is 4.
}

  \item{design}{
An optional formula for the DESeq2 model. If NULL, a default \code{"~batch + group"} design is used.
}

  \item{Ctrl_sample}{
String specifying the label of the control sample type used in DESeq2 normalisation. Default is \code{"normal"}.
}

  \item{rnaSeq}{
String specifying the type of RNA-seq data: \code{"bulk"} (default, uses DESeq2) or \code{"sc"} (single-cell, uses scuttle normalisation).
}

  \item{usePCA}{
Logical. Whether to reduce dimensionality with PCA before computing distance-based metrics. Default is TRUE.
}

  \item{nPCs}{
The number of principal components to use when \code{usePCA = TRUE}. Default is 50.
}

  \item{var.model}{
String specifying the count resampling model for \code{"resampling"} or \code{"resampleSamples"} variance estimation. Either \code{"multinom"} (default) or \code{"nbinom"}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is the main high-level wrapper function of the BatchMetrics package. The function first normalises the data (using DESeq2 VST for bulk RNA-seq or scuttle log-normalisation for single-cell), then for each requested variance estimation method it calls either \code{\link{QC_bootstrap}} or \code{\link{QC_resample}}. The result is a named list of quality control metric results for each variance estimation method.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A named list where each element corresponds to a requested variance estimation method (\code{var_measure}). Each element is the output of \code{\link{QC_bootstrap}} or \code{\link{QC_resample}}.
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
\code{\link{QC_bootstrap}}, \code{\link{QC_resample}}
}
\examples{
\dontrun{
#Load example count data provided with the package
data(counts)
data(batch)
data(group)

#Run the full QC wrapper with bootstrap variance estimation
results = QC_wrapper(counts, batch, group,
                     metrics = c("F-score", "Davies-Bouldin"),
                     iters = 10,
                     var_measure = "bootstrap",
                     cores = 2,
                     rnaSeq = "bulk")
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
