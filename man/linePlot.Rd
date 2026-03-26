\name{linePlot}
\alias{linePlot}
\alias{linePlot2}
\alias{plotResults}
\alias{plotMetrics}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Line plots for batch correction metric results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Functions for visualising quality control metric results from bootstrap or resampling experiments. \code{linePlot} creates a simple line plot with error bars. \code{linePlot2} creates a more detailed line plot with quantile-based ribbons for multiple signals. \code{plotResults} and \code{plotMetrics} are convenience wrappers for quickly plotting QC bootstrap results.
}
\usage{
linePlot(dat, title = "title")

linePlot2(result_data, filename = "linePlot.pdf", fileType = "pdf",
          metrics = NULL, order = "leftToRight", standardize = TRUE,
          quantiles = c(0.25, 0.5, 0.75),
          signals = c("bio", "batch", "ratio"))

plotResults(res_data,
            metrics = c("F-scores", "Davies-Bouldin index",
                        "Chi-square test", "kBET rejection rate"))

plotMetrics(object, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{x} here~~
A matrix or data frame of metric results, where columns alternate between means and standard deviations. Used by \code{linePlot}.
}

  \item{title}{
A character string used as the plot title. Default is \code{"title"} (used by \code{linePlot}).
}

  \item{result_data}{
A list of metric results (output of \code{\link{QC_bootstrap}} or similar), or a list of lists structured by signal type (used by \code{linePlot2} and \code{plotResults}).
}

  \item{filename}{
The path for the output file (used by \code{linePlot2} and \code{plotMetrics}).
}

  \item{fileType}{
String specifying the output file type: \code{"pdf"} (default) or \code{"png"} (used by \code{linePlot2}).
}

  \item{metrics}{
A character vector of metric names to plot. If NULL, all available metrics are plotted (used by \code{linePlot2} and \code{plotResults}). Note: for \code{plotResults}, these are display labels for plot titles (e.g., \code{"F-scores"}, \code{"Davies-Bouldin index"}), not the computation keys used by \code{\link{evalBatchEffect}}.
}

  \item{order}{
String specifying the order of the x-axis: \code{"leftToRight"} (default, largest correction level on the right) or \code{"rightToLeft"} (used by \code{linePlot2}).
}

  \item{standardize}{
Logical. Whether to standardise metric values to the [0, 1] range before plotting. Default is TRUE (used by \code{linePlot2}).
}

  \item{quantiles}{
A numeric vector of quantile levels to display as lines. Default is \code{c(0.25, 0.5, 0.75)} (used by \code{linePlot2}).
}

  \item{signals}{
A character vector specifying which signal types to plot (\code{"bio"}, \code{"batch"}, \code{"ratio"}). Default is all three (used by \code{linePlot2}).
}

  \item{res_data}{
A list of metric results with a \code{$means} component (used by \code{plotResults}).
}

  \item{object}{
An object returned by \code{\link{QC_bootstrap}} (used by \code{plotMetrics}).
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{linePlot} expects a matrix where odd-numbered columns contain metric means and even-numbered columns contain standard deviations. It produces a line plot with error bars.

\code{linePlot2} produces a more publication-ready plot showing quantile ranges across bootstrap or resampling iterations, saved to a PDF or PNG file.

\code{plotResults} is a convenience wrapper that calls \code{linePlot} for each metric.

\code{plotMetrics} normalises metric values and then calls \code{plotResults} to create a PNG file.
}
\value{
%%  ~Describe the value returned
\code{linePlot}: A \code{ggplot2} plot object.

\code{linePlot2}: A named list of \code{ggplot2} plot objects (one per metric), and saves the combined plot to the specified file.

\code{plotResults}: A list of \code{ggplot2} plot objects.

\code{plotMetrics}: None (saves a PNG to the path specified by \code{filename}).
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
\code{\link{QC_bootstrap}}, \code{\link{multiplot}}
}
\examples{
\dontrun{
#linePlot example: create a matrix of means and SDs
dat = matrix(rnorm(6 * 2), nrow = 3, ncol = 6)
p = linePlot(dat, title = "My metric")
print(p)
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
