\name{diluteSeries}
\alias{diluteSeries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Dilute Series
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Creates a list of datasets with controlled percentage of batch effects corrected for each dataset. Useful in evaluation of batch correction performance metrics.
}
\usage{
diluteSeries(Data,
batch1,
batch2 = NULL,
corrections1,
corrections2 = NULL,
perc = seq(0,1, by = 0.1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{x} here~~
The original dataset before batch correction
}

\item{batch1}{
The vector of labels of the first batch effect to be corrected.
}

\item{batch2}{
The vector of labels of the second batch effect to be corrected (optional, default is NULL).
}

\item{corrections1}{
The batch correction parameters for the first batch effect
}

\item{corrections2}{
The batch correction parameters for the second batch effect (optional, default is NULL).
}

\item{perc}{
The vector of percentages of batch correction to be added to each dataset in the list. Default is 0, 10\%,20\%, ..., 100\%.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A list of datasets with the chosen percentages of batch correction added
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

#Subset of data the calculate batch correction on
Subset = sample(1:100, size = 10, replace = FALSE)

#Perform batch correction
CorrData = BatchCorr_WithSubset(Data, batches, Subset, p = 1, SampleTypes=NULL,
                                 CorrMethod='combat', mean.only=FALSE)

#Extract the correction parameters
corrections = GetMeanCorrections(Data, CorrData, batches, DoSD=TRUE)

#Create a dilute series on those corrections
series = diluteSeries(Data, batch1 = batches, corrections1 = corrections)
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
