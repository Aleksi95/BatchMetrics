\name{ClusterCenter}
\alias{ClusterCenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Cluster centering
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A batch correction method that mean centers the batch clusters with additional option of variance correction.
}
\usage{
ClusterCenter(data, batch, mean.only = FALSE, p = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
The data matrix to be batch-corrected.
}

\item{batch}{
A vector of batch labels.
}

\item{mean.only}{
Logical whether to correct for variance also. Default is FALSE.
}

\item{p}{
The percentage of batch correction to do. Default is 1.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
The data matrix corrected for batch effects
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Developed by Aleksi Laiho (aleksi.laiho@helsinki.fi)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SamplesCenter}}
}
\examples{
#Generate data matrix and batches
Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
batches = sample(1:2, size = 100, replace = TRUE)

#Perform batch correction
ClusterCenter(Data, batches, mean.only = FALSE, p = 1)


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
