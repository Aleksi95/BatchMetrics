\name{getRanks}
\alias{getRanks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Average ranks of performance.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Returns the average ranks of performance of each method for all the metrics used for an object returned by \code{\link{QC_bootstrap}}, for the biological and batch categories, as well as their ratios.
}
\usage{
getRanks(metrics_list, type = c("bio", "batch", "ratio"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{metrics_list}{
%%     ~~Describe \code{x} here~~
An object returned by \code{\link{QC_bootstrap}}.
}

\item{type}{
The category of which the rank of metrics performance is calculated ("bio", "batch" and "ratio")
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
A list of data frames describing the performance rank for each method used
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

#calculate the ranks of best performing metrics
getRanks(qcmetrics, type = c("bio", "batch", "ratio"))
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
