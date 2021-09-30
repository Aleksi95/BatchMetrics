\name{BaLOOCV}
\alias{BaLOOCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Leave one out cross validation for batch correction
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Performs leave one out cross validation on batch correction methods and constructs batch-corrected dataset and gives out estimates for model performance. True sizes of batch effects may be given if known, in which case mean squared errors can be computed. In other cases, mean distance between batches is computed as a metric of batch correction performance.
}
\usage{
BaLOOCV(dataIn,
batch1,
batch2 = NULL,
sample_types = NULL,
subsetInd1 = NULL,
subsetInd2 = NULL,
method = "combat",
mse = FALSE,
true_params = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataIn}{
%%     ~~Describe \code{x} here~~
The gene expression data matrix for which batch correction should be performed in. Rows are genes, columns are samples.
}

\item{batch1}{
The vector of labels for the first batch effect.
}

\item{batch2}{
The vector of labels for the second batch effect (optional, default is NULL).
}

\item{sample_types}{
The vector of labels for the sample conditions of interest to be accounted for in batch correction (optional, default is NULL, needed for "harman", "subsetharman" and "limma")
}

\item{subsetInd1}{
The vector of indices of the subset of samples to be used for batch correction on the first batch effect. Needed for "subsetcombat" and "subsetharman", default is NULL.
}

\item{subsetInd2}{
The vector of indices of the subset of samples to be used for batch correction on the second batch effect. Needed for "subsetcombat" and "subsetharman", default is NULL.
}

\item{method}{
String denoting the batch correction method to be used in the cross-validation. One of "combat" (standard ComBat for a single batch), "combat2x" (ComBat for two batch effects), "samplescenter" (\code{\link{SamplesCenter}}), "clustercenter" (\code{\link{ClusterCenter}}), "harman", (Harman, Oytam et al 2016), "subsetcombat" (ComBat defined on a subset of data), "subsetharman" (harman defined on a subset) or "limma" (linear model from the "limma"-package).
}

\item{mse}{
Logical whether to compute the mean squared error. Only computed if the true batch effect size parameters are known. Default is FALSE
}

\item{true_params}{
The true batch effect parameters used for calculating the mean squared error. Default is NULL. Should be a list with the components "mu" and "delta".
}



}

\details{
%%  ~~ If necessary, more details than the description above ~~
The cross-validation aims to counter "over-correction" of methods such as ComBat, with the leave-one-out method for samples. The methods are applied on all but one sample, and the batch correction parameters are applied on the one sample left out. The performance of the batch correction method is evaluated as the mean distance of the sample to other batches. If the true size of batch effects are known in a randomized dataset for example, mean squared error of the calculated batch correction parameters is calculated.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{data}{The batch corrected dataset with each sample corrected separately}
%%  \item{comp2 }{Description of 'comp2'}

\item{dists}{The mean distance to other batches for each sample}
%% ...

\item{mu_mse}{The mean squared error for the batch correction shift parameters}

\item{delta_mse}{The mean squared error for the batch correction scale parameters}
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
#Generate a data matrix, batches and sample groups

Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
batches = sample(1:2, size = 100, replace = TRUE)
sample.groups = sample(1:4, size = 100, replace = TRUE)

#Perform the cross-validation
CV = BaLOOCV(dataIn = Data,
batch1 = batches,
batch2 = NULL,
sample_types = sample.groups,
subsetInd1 = NULL,
subsetInd2 = NULL,
method = "combat",
mse = FALSE,
true_params = NULL)

##for a controlled amount of batch effetcs
library('invgamma')


#Generate batch effects
mus = list()
deltas = list()
for(i in 1:length(unique(batches))){
  batchI = which(batches == unique(batches)[i])
  mu = rnorm(nrow(Data))
  delta = rinvgamma(nrow(Data), 1)
  mus[[i]] = mu
  rowsds = apply(Data[,batchI], 1, sd)
  deltas[[i]] = cbind(rowsds, delta)
  Data[,batchI] = delta*(Data[,batchI] + mu)
}

#Perform cross validation
CV = BaLOOCV(dataIn = Data,
batch1 = batches,
batch2 = NULL,
sample_types = sample.groups,
subsetInd1 = NULL,
subsetInd2 = NULL,
method = "combat",
mse = TRUE,
true_params = list(mu = mus, delta = deltas))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
\concept{ ~cross validation }
% \concept{ ~cpt2 }
% Use only one concept per line.
