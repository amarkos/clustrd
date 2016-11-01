\name{clusval}
\alias{clusval}

\title{
Distance-based statistics for cluster quality assessment.
}
\description{
This function computes two distance-based statistics (average silhouette widths and the Calinski-Harabasz index), which can be used for cluster quality assessment and decision about the number of clusters and dimensions.
}
\usage{
clusval(x, dst = "full")
}

\arguments{
\item{x}{An object of class \code{cluspca} or \code{clusmca}. For \code{cluspca} the distance measure used is the Euclidean distance, and for \code{clusmca} is Gower's distance}
 \item{dst}{Specifies the data used to compute the distances between objects. Options are full for the original data (after possible scaling) and low for the object scores in the low-dimensional space \code{(default = "full")}}
}

\value{
\item{ch}{Calinski-Harabasz index}
\item{asw}{Average silhouette width}
}
\seealso{
\code{\link{tune_cluspca}}, \code{\link{tune_clusmca}}
}
\examples{
data(USArrests, package = "datasets")
outCDR = cluspca(USArrests, 3, 2, alpha = 0.6, rotation = "varimax")
clusval(outCDR, dst = "full")$asw
}

