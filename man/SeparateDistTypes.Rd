\name{SeparateDistTypes}
\alias{SeparateDistTypes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Separate distance types
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Creates an adjacency matrix to separate between cluster and within cluster distances given a sample type vector
}
\usage{
SeparateDistTypes(SampleTypes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SampleTypes}{
%%     ~~Describe \code{x} here~~
  A vector of cluster labels (can be a factor, character or numeric)
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
  \item{type_matrix}{A matrix denoting distance types. Within cluster distances are denoted by numerical values greater than zero, within cluster distances by 0, with -1 on the diagonal}
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

#Generate sample type vector
SampleTypes = sample(1:2, size = 10, replace = TRUE)

#Get the separated distance type matrix
SeparateDistTypes(SampleTypes)
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
