\name{BatchCorr_WithSubset}
\alias{BatchCorr_WithSubset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Batch Correction on a subset of data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A wrapper for batch correction methods that calculate batch correction parameters based on a subset of data and applies them to.
}
\usage{
BatchCorr_WithSubset(DataIn, batch, SubsetIndex, p = 1, SampleTypes=NULL,
                                 CorrMethod='combat', mean.only=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataIn}{
%%     ~~Describe \code{x} here~~
The gene expression data to be batch corrected. Rows are genes and columns are samples.
}

\item{batch}{
The vector of batch labels.
}

\item{SubsetIndex}{
The indices of the columns to calculate the batch correction parameters on.
}

\item{p}{
The percentage of the batch correction parameters to be applied on the data (for \code{\link{diluteSeries}} for example). Default is 1 (100 \%).
}

\item{SampleTypes}{
The sample type vector to include in the batch correction analysis. Default is NULL.
}

\item{CorrMethod}{
The method of batch correction to calculate the batch correction parameters by ("combat", "harman", "samplecenter" or "clustercenter")
}

\item{mean.only}{
Logical whether to compute only the batch correction shift parameter or both shift and scale parameters. Default is FALSE.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function aims to counter the over-correction problem of batch effect correction. In particular this function can be utilised to calculate batch correction parameters based on technical replicate samples to counter for over-correction of biological variance that might distort downstream analyses. Utilises four different options for batch-correction methods chosen by the user, ComBat (Johnson et al, 2007), Harman (Oytam et al, 2016) \code{\link{SamplesCenter}} and \code{\link{ClusterCenter}}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
The batch-corrected data matrix
%% ...
}
\references{
%% ~put references to the literature/web site here ~
W. Evan Johnson, Cheng Li, Ariel Rabinovic, Adjusting batch effects in microarray expression data using empirical Bayes methods, Biostatistics, Volume 8, Issue 1, January 2007, Pages 118â127, https://doi.org/10.1093/biostatistics/kxj037

Oytam, Yalchin & Sobhanmanesh, Fariborz & Duesing, Konsta & Bowden, Joshua & Osmond-McLeod, Megan & Ross, Jason. (2016). Risk-conscious correction of batch effects: Maximising information extraction from high-throughput genomic datasets. BMC Bioinformatics. 17. 10.1186/s12859-016-1212-5.

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
\code{\link{SamplesCenter}} and \code{\link{ClusterCenter}}.
}
\examples{
#Generate a dataset and batches

Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
batches = sample(1:2, size = 100, replace = TRUE)

#Subset to calculate batch correction on
Subset = sample(1:100, size = 10, replace = FALSE)


#Perform batch correction
CorrData = BatchCorr_WithSubset(Data, batches, Subset, p = 1, SampleTypes=NULL,
                                 CorrMethod='combat', mean.only=FALSE)
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
