\name{plot.clusmca}
\alias{plot.clusmca}
\title{
Plotting function for \code{clusmca()} output.
}
\description{
Plotting function that creates a \code{ggplot2} based map of the object scores and a scatter plot of both the attribute scores and the centroids.
}
\usage{
\method{plot}{clusmca}(x, dims = c(1,2), disp = TRUE, cludesc = FALSE, what = c(TRUE,TRUE),
attlabs = NULL, binary = FALSE, ...)
}

\arguments{
  \item{x}{Object returned by \code{clusmca()}}
  \item{dims}{Numerical vector of length 2 indicating the dimensions to plot on horizontal and vertical axes respectively; default is first dimension horizontal and
second dimension vertical}
  \item{disp}{A logical value indicating whether the plots are shown in the \code{R} window or saved as PDF files in the working directory \code{(default = TRUE)}}
  \item{what}{Vector of two logical values specifying the contents of the plots. First entry indicates whether a scatterplot of the objects is displayed in principal coordinates. Second entry indicates whether a scatterplot of the attribute categories is displayed in principal coordinates. The default is \code{c(TRUE, TRUE)} and the resultant plot is a biplot of both objects and attribute categories with gamma-based scaling (see van de Velden et al. (2016))}
  \item{cludesc}{A logical value indicating whether a series of barplots is produced showing the largest (in absolute value) standardized residuals per attribute for each cluster \code{(default = FALSE)}}
  \item{attlabs}{Vector of attribute labels; if not provided, default labeling is applied}
  \item{binary}{Vector of attribute labels; if not provided, default labeling is applied}
  \item{...}{Further arguments to be transferred to \code{clusmca()}}
}

\references{
Hwang, H., Dillon, W. R. and Takane, Y. (2006). An extension of multiple correspondence analysis for identifying heterogenous subgroups of respondents, \emph{Psychometrika}, 71, 161-171.

Iodice D' Enza, A. and Palumbo, F. (2013). Iterative factor clustering of binary data. \emph{Computational Statistics}, \emph{28}(2), 789-807.

Van de Velden M., Iodice D' Enza, A. and Palumbo, F. (2016).  Cluster correspondence analysis.\emph{Psychometrika} (in press) DOI: 10.1007/s11336-016-9514-0
}
\seealso{
\link{plot.cluspca}
}
\examples{
data("hsq")
outclusMCA = clusmca(hsq[,1:8], 3, 2, method = "iFCB")
plot(outclusMCA, cludesc = TRUE)
}

