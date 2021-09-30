\name{SamplesCenter}
\alias{SamplesCenter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Sample Group Centering.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A batch correction method that aims to correct batch effect in chosen subsets/sample groups separately by centering the means in each group.
}
\usage{
SamplesCenter(data, batch, samples, p = 1)
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

\item{samples}{
A vector of labels of the sample groups of interest
}



\item{p}{
The percentage of batch correction to do.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
The batch correction should be performed for large enough sample groups (for too small groups there might be no batch effect or too little data). The method also does not correct for variance.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
The data matrix corrected for batch effects
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#Generate data matrix, batches and sample groups

Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
batches = sample(1:2, size = 100, replace = TRUE)
sample.groups = sample(1:4, size = 100, replace = TRUE)

#Perform batch correction
SamplesCenter(Data, batches, sample.groups, p = 1)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
\concept{ ~Batch correction }
\concept{ ~centering }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
