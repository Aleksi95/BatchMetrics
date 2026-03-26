\name{multiplot}
\alias{multiplot}
\alias{extract_legend}
\alias{std}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility functions for plot layout
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{multiplot} arranges multiple \code{ggplot2} plots in a grid layout. \code{extract_legend} extracts the legend from a \code{ggplot2} plot as a separate graphical object. \code{std} standardises rows of a matrix to the [0, 1] range.
}
\usage{
multiplot(..., plotlist = NULL, file, cols = 1, layout = NULL)

extract_legend(my_ggp)

std(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{x} here~~
\code{ggplot2} plot objects to arrange (used by \code{multiplot}).
}

  \item{plotlist}{
An optional list of \code{ggplot2} plot objects (used by \code{multiplot}).
}

  \item{file}{
Not currently used; retained for API compatibility (used by \code{multiplot}).
}

  \item{cols}{
The number of columns in the plot grid. Default is 1 (used by \code{multiplot}).
}

  \item{layout}{
An optional matrix specifying the layout of the plots. If NULL, the layout is determined automatically from \code{cols} (used by \code{multiplot}).
}

  \item{my_ggp}{
A \code{ggplot2} plot object from which to extract the legend (used by \code{extract_legend}).
}

  \item{d}{
A numeric matrix whose rows are to be standardised to the [0, 1] range (used by \code{std}).
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{multiplot} uses the \code{grid} package to arrange plots in a custom layout. If only one plot is provided it is printed directly, otherwise the grid layout is used.

\code{extract_legend} uses \code{ggplot_gtable} and \code{ggplot_build} to extract the guide box from a \code{ggplot2} plot, which can then be placed separately in a layout.

\code{std} applies min-max normalisation row-wise to a matrix so that each row ranges from 0 to 1.
}
\value{
%%  ~Describe the value returned
\code{multiplot}: None (prints the plots to the current graphics device).

\code{extract_legend}: A \code{grob} object containing the extracted legend.

\code{std}: A matrix with rows standardised to the [0, 1] range.
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
\code{\link{linePlot}}, \code{\link{linePlot2}}
}
\examples{
library(ggplot2)

p1 = ggplot(data.frame(x = 1:10, y = rnorm(10)), aes(x, y)) + geom_line()
p2 = ggplot(data.frame(x = 1:10, y = rnorm(10)), aes(x, y)) + geom_line()

\dontrun{
#Arrange two plots in two columns
multiplot(p1, p2, cols = 2)

#Extract the legend from a plot
p_legend = p1 + ggplot2::aes(color = factor(x))
leg = extract_legend(p_legend)
}

#Standardise a matrix row-wise
m = matrix(c(1, 3, 2, 8, 4, 10), nrow = 2)
std(m)
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
