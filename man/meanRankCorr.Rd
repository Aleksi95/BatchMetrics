\name{meanRankCorr}
\alias{meanRankCorr}
\alias{rankCorrelations}
\alias{rankCorrs2}
\alias{Ftest2}
\alias{writeLatexTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Metric evaluation utility functions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A collection of utility functions for evaluating and summarising the performance of batch correction quality metrics. \code{meanRankCorr} computes the mean Spearman rank correlation between metric values and a rank vector. \code{rankCorrelations} computes rank correlations for multiple metrics across dilution levels. \code{rankCorrs2} computes Spearman rank correlations for multiple metrics and signal types simultaneously. \code{Ftest2} tests whether metric values differ significantly across dilution levels. \code{writeLatexTable} writes a matrix or data frame to a LaTeX tabular environment.
}
\usage{
meanRankCorr(resdata, rank = 1:11)

rankCorrelations(resdata, type = "spearman", score = "ratio",
                 levels = seq(0, 1, by = 0.1), metrics = NULL)

rankCorrs2(resdata, metrics = NULL, types = NULL)

Ftest2(data, test = "AOV", alternative = "increasing")

writeLatexTable(table, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resdata}{
%%     ~~Describe \code{x} here~~
A list of metric results, typically from \code{\link{QC_bootstrap}} or \code{\link{QC_resample}}, structured by signal type and metric name. Used by \code{meanRankCorr}, \code{rankCorrelations} and \code{rankCorrs2}.
}

  \item{rank}{
A numeric vector of expected ranks to correlate against. Default is \code{1:11} (corresponding to the 11 dilution levels from 0 to 1 by 0.1). Used by \code{meanRankCorr}.
}

  \item{type}{
String specifying the correlation type: \code{"spearman"} (default) or \code{"pearson"}. Used by \code{rankCorrelations}.
}

  \item{score}{
String specifying which signal type to use for correlation: \code{"bio"}, \code{"batch"} or \code{"ratio"} (default). Used by \code{rankCorrelations}.
}

  \item{levels}{
A numeric vector of the dilution levels used in the experiment. Default is \code{seq(0, 1, by = 0.1)}. Used by \code{rankCorrelations}.
}

  \item{metrics}{
A character vector of metric names to include. If NULL, all available metrics are used. Used by \code{rankCorrelations} and \code{rankCorrs2}.
}

  \item{types}{
A character vector of signal type names to include. If NULL, all available signal types are used. Used by \code{rankCorrs2}.
}

  \item{data}{
A matrix or data frame where rows are dilution levels and columns are bootstrap/resampling iterations. Used by \code{Ftest2}.
}

  \item{test}{
String specifying the statistical test to use: \code{"AOV"} (ANOVA, default) or \code{"JT"} (Jonckheere-Terpstra trend test). Used by \code{Ftest2}.
}

  \item{alternative}{
String specifying the alternative hypothesis direction for the Jonckheere-Terpstra test: \code{"increasing"} (default), \code{"decreasing"} or \code{"two.sided"}. Used by \code{Ftest2}.
}

  \item{table}{
A matrix or data frame to write as a LaTeX table. Used by \code{writeLatexTable}.
}

  \item{filename}{
The path to the output file (used by \code{writeLatexTable}).
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{meanRankCorr} computes the Spearman correlation between each column of the transposed metric data and the expected rank vector, then returns the mean of these correlations.

\code{rankCorrelations} pools all bootstrap/resampling values for each metric across dilution levels and computes a Spearman or Pearson correlation with the corresponding dilution level values.

\code{rankCorrs2} computes Spearman correlations for all metric-by-signal-type combinations, with batch-type correlations negated so that higher values consistently mean better batch correction performance.

\code{Ftest2} tests whether metric values vary systematically across dilution levels using a one-way ANOVA (F-test) or the Jonckheere-Terpstra ordered alternative test.

\code{writeLatexTable} writes a table to a plain-text file formatted as a LaTeX \code{tabular} environment with horizontal rules.
}
\value{
%%  ~Describe the value returned
\code{meanRankCorr}: A single numeric value (mean Spearman correlation).

\code{rankCorrelations}: A named list with components:
\itemize{
  \item \code{estimates}: A named numeric vector of correlation estimates for each metric.
  \item \code{p.values}: A named numeric vector of p-values for each metric.
}

\code{rankCorrs2}: A numeric matrix with metrics as columns and signal types as rows, containing Spearman correlation coefficients.

\code{Ftest2}: For \code{test = "AOV"}, a single numeric F-statistic. For \code{test = "JT"}, a named numeric vector with components \code{JT} (the test statistic) and \code{pvalue}.

\code{writeLatexTable}: None (writes the LaTeX table to \code{filename}).
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
\code{\link{QC_bootstrap}}
}
\examples{
#Create a simple matrix (rows = dilution levels, cols = iterations)
m = matrix(c(seq(1, 11), seq(11, 1, by = -1),
             sample(1:11, 11)), nrow = 11, ncol = 3)

#F-test across groups
Ftest2(m, test = "AOV")

#Write a small table to a LaTeX file
tab = data.frame(Method = c("ComBat", "Harman"),
                 Score = c(0.85, 0.72))
\dontrun{
writeLatexTable(tab, filename = tempfile(fileext = ".tex"))
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
