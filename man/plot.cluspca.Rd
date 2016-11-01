\name{plot.cluspca}
\alias{plot.cluspca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting function for \code{cluspca()} output.
}
\description{
Plotting function that creates a \code{ggplot2} scatterplot of the objects, a correlation circle of the variables or a biplot of both objects and variables.}

\usage{
\method{plot}{cluspca}(x, dims = c(1, 2), disp = TRUE, cludesc = FALSE, what = c(TRUE,TRUE), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object returned by \code{cluspca()}}
  \item{dims}{Numerical vector of length 2 indicating the dimensions to plot on horizontal and vertical axes respectively; default is first dimension horizontal and
second dimension vertical}
  \item{disp}{A logical value indicating whether the plots are shown in the \code{R} window or saved as PDF files in the working directory \code{(default = TRUE)}}
  \item{what}{Vector of two logical values specifying the contents of the plots. First entry indicates whether a scatterplot of the objects is displayed and the second entry whether a correlation circle of the variables is displayed. The default is \code{c(TRUE, TRUE)} and the resultant plot is a biplot of both objects and variables}
  \item{cludesc}{A logical value indicating if a parallel plot showing cluster means is produced (default = \code{FALSE)}}
  \item{...}{Further arguments to be transferred to \code{cluspca()}}
}
\references{
De Soete, G. and Carroll, J. D. (1994). K-means clustering in a low-dimensional Euclidean space. In Diday E. et al. (Eds.), \emph{New Approaches in Classification and Data Analysis}, Heidelberg: Springer, 212-219.

Vichi, M. and Kiers, H.A.L. (2001). Factorial K-means analysis for two-way data. \emph{Computational Statistics and Data Analysis}, 37, 49-64.
}
\seealso{
\code{\link{plot.clusmca}}
}
\examples{
data("iris", package = "datasets")
outclusPCA = cluspca(iris[,-5], 3, 2, alpha = 0.3, rotation = "varimax")
table(outclusPCA$cluID,iris[,5])
plot(outclusPCA, cludesc = TRUE)
}

