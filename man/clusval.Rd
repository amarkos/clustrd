\name{clusval}
\alias{clusval}

\title{
Distance-based statistics for cluster quality assessment.
}
\description{
This function computes two distance-based statistics (average silhouette width and the Calinski-Harabasz index), which can be used for cluster quality assessment and decision about the number of clusters and dimensions.
}
\usage{
clusval(x, dst = "full")
}

\arguments{
\item{x}{An object of class \code{cluspca} or \code{clusmca}. For \code{cluspca} the distance measure used is the Euclidean distance, and for \code{clusmca} is Gower's distance}
 \item{dst}{Specifies the data used to compute the distances between objects. Options are full for the original data (after possible scaling) and low for the object scores in the low-dimensional space (default = \code{"full")}}
}

\value{
\item{ch}{Calinski-Harabasz index}
\item{asw}{Average silhouette width}
}
\references{
Calinski, R.B., and Harabasz, J., (1974). A dendrite method for cluster analysis. \emph{Communications in Statistics}, 3, 1-27.

Kaufman, L., and Rousseeuw, P.J., (1990). \emph{Finding Groups in Data: An Introduction to Cluster Analysis}. Wiley, New York.
}

\seealso{
\code{\link{tuneclus}}, \code{\link{cluspca}}, \code{\link{clusmca}}
}
\examples{
data(USArrests, package = "datasets")
#Reduced K-means on the USArrests dataset 
#4 clusters in 3 dimensions after 100 random starts
outclusPCA = cluspca(USArrests, 4, 3, rotation = "varimax")
#Cluster quality assessment in the full dimensional space
clusval(outclusPCA, dst = "full")
}

