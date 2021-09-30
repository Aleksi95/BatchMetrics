\name{ApplyBatchCorrections}
\alias{ApplyBatchCorrections}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Apply batch correction parameters
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Applies batch corrections parameters computed by \code{\link{GetMeanCorrections}} to a dataset.
}
\usage{
ApplyBatchCorrections(DataIn, Batches, corrections, p=1, mean.only = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataIn}{
%%     ~~Describe \code{x} here~~
The data matrix before batch corrections to apply the batch correction parameters on to.
}

\item{Batches}{
The vector of batch labels
}

\item{corrections}{
The list of batch correction parameters produced by \code{\link{GetMeanCorrections}}.
}

\item{p}{
The percentage of batch correction parameters to apply on the data (for example for \code{\link{diluteSeries}})
}

\item{mean.only}{
Logical whether to only apply mean corrections.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
The batch corrected data matrix
%%  \item{comp2 }{Description of 'comp2'}
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
#Data matrix
Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
batches = sample(1:2, size = 100, replace = TRUE)

#Subset to calculate batch correction on
Subset = sample(1:100, size = 10, replace = FALSE)

#Perform batch correction
CorrData = BatchCorr_WithSubset(Data, batches, Subset, p = 1, SampleTypes=NULL,
                                 CorrMethod='combat', mean.only=FALSE)

#Extract batch correction parameters
corrections = GetMeanCorrections(Data, CorrData, batches, DoSD=TRUE)

#Apply batch correction parameters on the data
CorrectedData = ApplyBatchCorrections(Data, batches, corrections)
#Should be the same dataset as CorrData

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
