\name{FindMatchIndex}
\alias{FindMatchIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Find matcg index
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Internal utility function. Which cases of vector1 are found in vector2 and return an index vector of these positions.
}
\usage{
FindMatchIndex(vector1, vector2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vector1}{
%%     ~~Describe \code{x} here~~
First vector
}

\item{vector2}{
Second vector
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
%% ...
Index vector of matching positions.
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

\dontrun{FindMatchIndex(vector1, vector2)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\keyword{ ~internal }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
