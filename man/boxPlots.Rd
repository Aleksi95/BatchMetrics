\name{boxPlots}
\alias{boxPlots}
\alias{boxPlot}
\alias{boxPlotPDF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Boxplot of bootstrap results
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Produces a boxplot of the quality control metric results for one or multiple datasets
}
\usage{
boxPlots(object)
boxPlotPDF(object, filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{x} here~~
An object returned by \code{\link{QC_bootstrap}}

}
\item{filename}{
A filename for the pdf returned by \code{\link{boxPlotPDF}}. Must end with ".pdf".
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
boxPlots prints the plots in an R environment, boxPlotPDF prints the boxplots as a PDF in the user's working directory.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
None
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

#Plot the results
\dontrun{boxPlots(qcmetrics)}

#Plot the results to PDF
\dontrun{boxPlotPDF(qcmetrics, "metrics.pdf")}
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
