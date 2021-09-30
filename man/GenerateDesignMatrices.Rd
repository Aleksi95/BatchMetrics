\name{GenerateDesignMatrices}
\alias{GenerateDesignMatrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Generate Design Matrices
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generates a design matrix of the sample conditions.
}
\usage{
GenerateDesignMatrices(SampleClassifications)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SampleClassifications}{
%%     ~~Describe \code{x} here~~
A vector of labels of the sample classifications of interest.
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
A design matrix of the desired sample conditions.
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

sample.groups = sample(1:4, size = 100, replace = TRUE)

GenerateDesignMatrices(sample.groups)
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
