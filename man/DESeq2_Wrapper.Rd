\name{DESeq2_Wrapper}
\alias{DESeq2_Wrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DESeq2 wrapper for normalisation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A convenience wrapper around the \code{DESeq2} workflow that constructs a \code{DESeqDataSet}, estimates size factors, and estimates dispersions from a raw count matrix.
}
\usage{
DESeq2_Wrapper(CountData, coldata, designCols, ctrl_sample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CountData}{
%%     ~~Describe \code{x} here~~
A raw RNA-seq count matrix with rows as genes and columns as samples.
}

  \item{coldata}{
A data frame of sample annotations. Column names must include the variables referenced in \code{designCols}.
}

  \item{designCols}{
A character vector of column names from \code{coldata} to include in the DESeq2 design formula (e.g., \code{c("batch", "group")}).
}

  \item{ctrl_sample}{
String specifying the control/reference sample type (currently not directly used in the formula construction but retained for compatibility).
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function builds a design formula by concatenating the \code{designCols} with \code{+}, creates a \code{DESeqDataSet}, estimates size factors (quietly) and estimates dispersions (quietly). The returned object can be passed to downstream DESeq2 functions such as \code{vst}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A \code{DESeqDataSet} object with size factors and dispersions estimated.
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, 550. https://doi.org/10.1186/s13059-014-0550-8
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
\code{\link{QC_wrapper}}, \code{\link{QC_resample}}
}
\examples{
\dontrun{
#Load example count data provided with the package
data(counts)
data(batch)
data(group)

coldata = data.frame(batch = batch, group = group)
rownames(coldata) = colnames(counts)

#Run DESeq2 wrapper
deseq_obj = DESeq2_Wrapper(counts, coldata,
                            designCols = c("batch", "group"),
                            ctrl_sample = "normal")
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
