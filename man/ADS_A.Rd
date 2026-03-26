\name{ADS_A}
\alias{ADS_A}
\alias{diluteSeries2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample-based dilute series for batch effect simulation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{ADS_A} removes a fraction \code{p} of the sample-group mean signal from the data, effectively diluting biological group differences. \code{diluteSeries2} is a wrapper that applies \code{ADS_A} across a vector of dilution levels to generate a list of datasets with varying degrees of biological group separation.
}
\usage{
ADS_A(data, samples, p = 1)

diluteSeries2(data, samples, perc = seq(0, 1, by = 0.1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
A gene expression data matrix with rows as features (genes) and columns as samples.
}

  \item{samples}{
A vector of sample group labels (biological groups) for each sample.
}

  \item{p}{
The proportion of the group mean to subtract. A value of 1 removes all between-group mean differences; a value of 0 leaves the data unchanged. Default is 1.
}

  \item{perc}{
A numeric vector of dilution levels passed to \code{ADS_A}. Default is \code{seq(0, 1, by = 0.1)}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{ADS_A} subtracts a fraction \code{p} of the within-group column means from all columns in each group. At \code{p = 1}, the within-group means are entirely removed, eliminating systematic between-group differences. At \code{p = 0}, the data is unchanged.

\code{diluteSeries2} applies \code{ADS_A} for each level in \code{perc} and returns the resulting list of data matrices. This is useful for simulating datasets with varying degrees of biological signal when evaluating the sensitivity of batch correction quality metrics.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
For \code{ADS_A}: A data matrix of the same dimensions as \code{data}, with the group mean signal diluted by the factor \code{p}.

For \code{diluteSeries2}: A list of data matrices, one for each dilution level in \code{perc}.
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
\code{\link{diluteSeries}}
}
\examples{
#Generate a gene expression data matrix (genes x samples)
Data = matrix(rnorm(200 * 20, mean = 0, sd = 1), nrow = 200, ncol = 20)

#Define biological group labels
samples = sample(c("A", "B"), 20, replace = TRUE)

#Dilute 50% of the group mean differences
diluted = ADS_A(Data, samples, p = 0.5)
dim(diluted)

#Generate a full dilute series
series = diluteSeries2(Data, samples, perc = seq(0, 1, by = 0.25))
length(series)
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
