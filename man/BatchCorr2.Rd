\name{BatchCorr2}
\alias{BatchCorr2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Batch correction for two batch effects
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A wrapper for two different batch correction methods to correct for two batch effects.
}
\usage{
BatchCorr2(dataIn,
            batch1,
            batch2,
            samples = NULL,
            subset1 = NULL,
            subset2 = NULL,
            p = 1,
            method = "combat")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataIn}{
%%     ~~Describe \code{x} here~~
The data matrix to correct batch effects for.
}

\item{batch1}{
The vector of labels of the first batch effect to be corrected.
}

\item{batch2}{
The vector of labels of the second batch effect to be corrected.
}

\item{samples}{
The vector of labels for the sample groups of interest to be accounted for in batch correction (optional, default is NULL, needed for Harman and samplescenter)
}

\item{subset1}{
The subset of technical replicates in the first batch to compute the batch correction parameters for (optional).
}

\item{subset2}{
The subset of technical replicates in the second batch to compute the batch correction parameters for (optional).
}

\item{p}{
The percentage of batch correction to do. Default is 1.
}

\item{method}{
String denoting the batch correction method to use ("combat", "harman", "clustercenter" or "samplescenter"). Note: if subset1 and subset2 is not null, only "combat" and "harman" are supported. "harman" and "samplescenter" also require sample type vector to be used.
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
The data matrix corrected for batch effects
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
#Generate data matrix, batches and sample groups

Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
batch1 = sample(1:2, size = 100, replace = TRUE)
batch2 = sample(1:2, size = 100, replace = TRUE)

sample.groups = sample(1:4, size = 100, replace = TRUE)

#Perform batch correction with 2 batch effects
corrData = BatchCorr2(Data,
batch1,
batch2,
samples = sample.groups,
subset1 = NULL,
subset2 = NULL,
p = 1,
method = "combat")
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
