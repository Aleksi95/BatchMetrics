\name{GetMeanCorrections}
\alias{GetMeanCorrections}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Extract batch corrections.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Extracts the average shift and scale parameters created in batch correction.
}
\usage{
GetMeanCorrections(OldData, NewData, batch, DoSD=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{OldData}{
%%     ~~Describe \code{x} here~~
The data matrix before batch correction.
}

\item{NewData}{
The data matrix after batch correction.
}

\item{batch}{
The vector of batch labels.
}

\item{DoSD}{
Logical whether to compute average scale parameter of batch correction along with shift parameters. Default is TRUE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
Extracting average batch correction parameters is utilised in functions such as \code{\link{BatchCorr_WithSubset}} and \code{\link{diluteSeries}}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{mean.vals}{The average batch correction shift parameters for each batch}
\item{sd.vals}{The average batch correction scale parameters for each batch}
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
\code{\link{BatchCorr_WithSubset}} and \code{\link{diluteSeries}}
}
\examples{
#Generate dataset and batches

Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
batches = sample(1:2, size = 100, replace = TRUE)

#Subset to calculate batch correction on
Subset = sample(1:100, size = 10, replace = FALSE)

#Perform batch correction
CorrData = BatchCorr_WithSubset(Data, batches, Subset, p = 1, SampleTypes=NULL,
                                 CorrMethod='combat', mean.only=FALSE)

#Extract batch correction parameters
corrections = GetMeanCorrections(Data, CorrData, batches, DoSD=TRUE)
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
