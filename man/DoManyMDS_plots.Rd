\name{DoManyMDS_plots}
\alias{DoManyMDS_plots}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Multiple multidimensional scaling plots
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Produces multiple multidimensional scaling plots (see \code{\link{plotMDS_andLegend}}), given the selected columns of the sample type table as legend.
}
\usage{
DoManyMDS_plots(data,
SampleTypeTable,
ColVector,
SelectScores = NULL,
toPDF = FALSE,
filename = NULL,
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
The given gene expression dataset, columns are samples, rows are genes.
}

\item{SampleTypeTable}{
The data table containing the information of the sample columns in the data, with the columns containing the labels for the desired groups.
}

\item{ColVector}{
The vector denoting the columns of SampleTypeTable which to select.
}

\item{SelectScores}{
Vector of scores the order of which the genes should be selected by (see \code{\link{plotMDS_andLegend}}).
}

\item{toPDF}{
Logical whether to print the plot as a PDF file.
}

\item{filename}{
The beginning string of the filename of the pdf. The final filename will be of the form "filename_MDSplot.pdf". Default is NULL
}

\item{...}{
Additional arguments
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
%% ...
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
\code{\link{plotMDS_andLegend}}
}
\examples{
library(limma)
#Generate dataset and batches

Data = matrix(rnorm(100*1000,mean=0,sd=1), ncol=100, nrow=1000)
batches = sample(1:2, size = 100, replace = TRUE)
sampletypes = sample(1:4, size = 100, replace = TRUE)

#Subset to calculate batch correction on
Subset = sample(1:100, size = 10, replace = FALSE)

#Perform batch correction
CorrData = BatchCorr_WithSubset(Data, batches, Subset, p = 1, SampleTypes=sampletypes,
                                 CorrMethod='combat', mean.only=FALSE)

#Plot results
DoManyMDS_plots(Data, SampleTypeTable = data.frame(batches, sampletypes),
ColVector = c(1,2), SelectScores = NULL, toPDF = FALSE, filename = NULL)
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
