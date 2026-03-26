\name{gPCA.batchdetect}
\alias{gPCA.batchdetect}
\alias{gPCA_percentage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Guided PCA for batch effect detection
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{gPCA.batchdetect} performs a guided principal component analysis (gPCA) to detect batch effects using a permutation test. \code{gPCA_percentage} is a convenience wrapper that returns the percentage of variation explained by batch relative to the first guided PC.
}
\usage{
gPCA.batchdetect(x, batch, filt = NULL, nperm = 1000, center = FALSE,
                 scaleY = FALSE, seed = NULL)

gPCA_percentage(data, sample_types, nperm = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
A numeric matrix with rows as samples and columns as features (genes). Used by \code{gPCA.batchdetect}.
}

  \item{batch}{
A numeric vector of batch labels (integer codes) for each sample. Used by \code{gPCA.batchdetect}.
}

  \item{filt}{
Optional integer specifying the number of top features (by standard deviation) to retain after filtering. Default is NULL (no filtering). Used by \code{gPCA.batchdetect}.
}

  \item{nperm}{
The number of permutations to perform during the permutation test. Default is 1000.
}

  \item{center}{
Logical. Indicates whether the data is already centered. If FALSE (default), the function centers the data using \code{scale(x, center = TRUE, scale = FALSE)} before PCA. If TRUE, the data is assumed to already be centered and is used as-is.
}

  \item{scaleY}{
Logical. Whether to scale the batch indicator matrix by the number of samples in each batch. Default is FALSE.
}

  \item{seed}{
Optional integer to set the random seed for reproducibility. Default is NULL.
}

  \item{data}{
A gene expression data matrix with rows as features (genes) and columns as samples. Used by \code{gPCA_percentage}.
}

  \item{sample_types}{
A vector of labels (factor, character or numeric) describing the batch or biological groups of each sample. Used by \code{gPCA_percentage}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
The gPCA method contrasts the first principal component of an unguided (standard) PCA with the first principal component of a guided PCA, where the PCA is guided by a batch indicator matrix. The test statistic delta is the ratio of the variance of the first guided PC to the variance of the first unguided PC. A permutation test is used to compute a p-value.

\code{gPCA_percentage} converts the sample type vector to numeric batch codes and returns the relative percentage of variance captured by batch in the guided PCA:
\deqn{(varPCg1 - varPCu1) / varPCg1}
}
\value{
For \code{gPCA.batchdetect}: A list with the following components:
\item{delta}{The observed test statistic (ratio of guided to unguided PC1 variances).}
\item{p.val}{The permutation p-value.}
\item{delta.p}{The vector of permuted delta values.}
\item{batch}{The batch vector used.}
\item{filt}{The number of features retained (or NULL if no filtering).}
\item{n}{Number of samples.}
\item{p}{Number of features.}
\item{b}{Number of batches.}
\item{PCg}{The guided PC scores.}
\item{PCu}{The unguided PC scores.}
\item{varPCu1}{Variance explained by the first unguided PC.}
\item{varPCg1}{Variance explained by the first guided PC.}
\item{nperm}{Number of permutations performed.}
\item{cumulative.var.u}{Cumulative variance explained by unguided PCs.}
\item{cumulative.var.g}{Cumulative variance explained by guided PCs.}

For \code{gPCA_percentage}: A single numeric value representing the percentage of variation in the first guided PC that is in excess of the first unguided PC.
}
\references{
%% ~put references to the literature/web site here ~
Reese, S. E., Archer, K. J., Therneau, T. M., Atkinson, E. J., Vachon, C. M., de Andrade, M., Kocher, J.-P. A., & Eckel-Passow, J. E. (2013). A new statistic for identifying batch effects in high-throughput genomic data that uses guided principal component analysis. Bioinformatics, 29(22), 2877-2883. https://doi.org/10.1093/bioinformatics/btt480
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
The permutation test in \code{gPCA.batchdetect} can be slow for large datasets with many permutations. Use \code{seed} for reproducibility and reduce \code{nperm} for faster exploration.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{EffectSizeEstimate}}
}
\examples{
#Generate a data matrix (samples x features, as required by gPCA.batchdetect)
Data = matrix(rnorm(20 * 200, mean = 0, sd = 1), nrow = 20, ncol = 200)

#Define integer batch labels
batch = sample(1:2, 20, replace = TRUE)

#Detect batch effects with gPCA (using few permutations for speed)
result = gPCA.batchdetect(Data, batch, nperm = 50, seed = 42)
result$delta
result$p.val

#Using the convenience wrapper (genes x samples format)
DataT = t(Data)
sample_types = sample(c("A", "B"), 20, replace = TRUE)
\dontrun{
pct = gPCA_percentage(DataT, sample_types, nperm = 50)
pct
}
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
