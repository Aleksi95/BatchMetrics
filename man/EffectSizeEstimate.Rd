\name{EffectSizeEstimate}
\alias{EffectSizeEstimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effect size estimate for batch or biological group separation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Estimates the effect size of batch or biological group separation by performing dimensionality reduction (MDS or PCA) and fitting a linear model.
}
\usage{
EffectSizeEstimate(data, sample_types, dim.reduct.method = "MDS")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
A gene expression data matrix with rows as features (genes) and columns as samples.
}

  \item{sample_types}{
A vector of labels (factor, character or numeric) describing the batch or biological group of each sample.
}

  \item{dim.reduct.method}{
String specifying the dimensionality reduction method to use. Either \code{"MDS"} (default, uses limma's \code{plotMDS}) or \code{"PCA"} (uses \code{prcomp}).
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function reduces the data to a one-dimensional summary by averaging the first two dimensions from the chosen dimensionality reduction method. A linear model is then fitted regressing this summary against the sample type labels. The function returns the estimated effect size as the absolute value of the regression coefficient, the p-value from the linear model, and the R-squared value.
}
\value{
%%  ~Describe the value returned
A named list with the following components:
\item{coef}{The absolute value of the regression coefficient (effect size estimate).}
\item{pval}{The p-value from the linear model test for the group effect.}
\item{R2}{The R-squared value of the linear model.}
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
\code{\link{gPCA.batchdetect}}, \code{\link{gPCA_percentage}}
}
\examples{
#Generate a gene expression data matrix (genes x samples)
Data = matrix(rnorm(1000 * 20, mean = 0, sd = 1), nrow = 1000, ncol = 20)

#Define sample group labels
sample_types = sample(c("group1", "group2"), 20, replace = TRUE)

#Estimate the effect size using MDS
result = EffectSizeEstimate(Data, sample_types, dim.reduct.method = "MDS")
result$coef
result$pval
result$R2
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
