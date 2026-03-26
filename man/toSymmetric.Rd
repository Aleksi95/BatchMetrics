\name{toSymmetric}
\alias{toSymmetric}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a triangular matrix to a symmetric matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Fills in the missing triangle of a triangular matrix to produce a symmetric matrix.
}
\usage{
toSymmetric(dist, lower = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{
%%     ~~Describe \code{x} here~~
A square numeric matrix that contains values in only one triangle (upper or lower). The diagonal can be zero or any value.
}

  \item{lower}{
Logical. If TRUE (default), the matrix contains data in the lower triangle and the upper triangle is filled in. If FALSE, the upper triangle is mirrored to the lower triangle.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
This is an internal utility function used to convert a lower or upper triangular distance matrix into a full symmetric matrix by mirroring the filled triangle.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A symmetric square matrix.
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
\code{\link{GetDistMatrix}}
}
\examples{
#Create a lower triangular matrix
m = matrix(0, nrow = 4, ncol = 4)
m[lower.tri(m)] = c(1, 2, 3, 4, 5, 6)

#Convert to a symmetric matrix
sym_m = toSymmetric(m, lower = TRUE)
sym_m
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
