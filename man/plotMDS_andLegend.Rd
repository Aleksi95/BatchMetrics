\name{plotMDS_andLegend}
\alias{plotMDS_andLegend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Multidimensional scaling plot with legend
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Plots the gene expression dataset with the given sample types and batch effects highlighted in the legend. Uses the \code{\link{limma}} package for the multidimensional scaling.
}
\usage{
plotMDS_andLegend(data, sample_types, LegendChoice="bottomleft",
                              SelectScores = NULL, MaxNum = 500, UseMax=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
The given gene expression dataset, columns are samples, rows are genes.
}

\item{sample_types}{
The vector of labels of the sample types of interest (batch effects or biological groups)
}

\item{LegendChoice}{
The string denoting the position of the legend, such as "bottomleft" or "bottomright" (see \code{\link{legend}}). Default is "bottomright".
}

\item{SelectScores}{
Vector of scores the order of which the genes should be selected by (in decreasing order, if UseMax = TRUE), such as the variance. Default is NULL. Has to be the same length as the number of rows in the data.
}

\item{MaxNum}{
The maximum number of genes to select. Default is 500.
}

\item{UseMax}{
Logical whether to select the genes in decreasing order given by the chosen score in SelectScores. Default is TRUE.
}

\item{...}{
Additional arguments.
}






}

\references{
%% ~put references to the literature/web site here ~
Cox, T.F.; Cox, M.A.A. (2001). Multidimensional Scaling. Chapman and Hall
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
\code{\link{DoManyMDS_plots}}
}
\examples{
#Generate dataset and batches

Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
batches = sample(1:2, size = 100, replace = TRUE)


#Subset to calculate batch correction on
Subset = sample(1:100, size = 10, replace = FALSE)

#Perform batch correction
CorrData = BatchCorr_WithSubset(Data, batches, Subset, p = 1, SampleTypes=NULL,
                                 CorrMethod='combat', mean.only=FALSE)

#Plot results
plotMDS_andLegend(Data, batches, LegendChoice="bottomleft",
                              SelectScores = NULL, MaxNum = 500, UseMax=TRUE, ...)
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
