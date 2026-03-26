\name{getCMS}
\alias{getCMS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cell Mixing Score (CMS)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Computes the mean Cell Mixing Score (CMS) from the CellMixS package as a metric for batch effect evaluation. Lower values indicate better mixing between batches.
}
\usage{
getCMS(data, batchvar, k = 50, doPCA = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
A gene expression data matrix with rows as features (genes) and columns as samples. If \code{doPCA = TRUE}, this should be a PCA-reduced matrix (rows as PCs, columns as samples).
}

  \item{batchvar}{
A vector of batch labels (factor, character or numeric) describing the batch or biological group of each sample.
}

  \item{k}{
The number of nearest neighbours to use in the CMS computation. Default is 50.
}

  \item{doPCA}{
Logical. Whether PCA has already been performed on the data (i.e., the data argument contains PCA scores). Default is FALSE.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
This function wraps the \code{cms} function from the \code{CellMixS} package. The CMS evaluates the batch mixing by computing a local test statistic for each cell, comparing the composition of its k-nearest neighbourhood to the global batch composition. The mean CMS across all cells is returned. The function internally creates a \code{SingleCellExperiment} object with the data and batch variable.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A single numeric value representing the mean CMS across all samples. Lower values indicate better batch mixing.
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Lütge, A., Zyprych-Walczak, J., Kunzmann Brydon, U., Hintermayr, M., Nowicka, M., Böhne, W., Robinson, M. D., & Soneson, C. (2021). CellMixS: quantifying and visualizing batch effects in single-cell RNA-seq data. Life Science Alliance, 4(6), e202001004. https://doi.org/10.26508/lsa.202001004
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
\code{\link{getkBET}}, \code{\link{getKNNprop}}
}
\examples{
#Generate a small expression data matrix (genes x samples)
Data = matrix(rnorm(50 * 20, mean = 0, sd = 1), nrow = 50, ncol = 20)

#Define batch labels
batchvar = sample(c("batch1", "batch2"), 20, replace = TRUE)

\dontrun{
#Compute the mean CMS
cms_score = getCMS(Data, batchvar, k = 10, doPCA = FALSE)
cms_score
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
