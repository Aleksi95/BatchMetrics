\name{removeZeroSD}
\alias{removeZeroSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Remove zero standard deviations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Removes standard zero valued standard deviations from correction parameter data
}
\usage{
removeZeroSD(corr_data, batches)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{corr_data}{
%%     ~~Describe \code{x} here~~
Correction parameter data produced by \code{\link{GetMeanCorrections}}.
}

\item{batches}{
Vector of batch labels.
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
Correction data with zero SDs removed.
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

\dontrun{
removeZeroSD()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{~Internal}
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
