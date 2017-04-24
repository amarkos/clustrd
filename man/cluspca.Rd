\name{cluspca}
\alias{cluspca}

\title{
Joint  dimension  reduction  and  clustering  of  continuous  data.
}
\description{
This function implements Factorial K-means (Vichi and Kiers, 2001) and Reduced K-means (De Soete and Carroll, 1994), as well as a compromise version of these two methods. The methods combine Principal Component Analysis for dimension reduction with K-means for clustering.
}

\usage{
cluspca(data, nclus, ndim, alpha = NULL, method = "RKM", 
center = TRUE, scale = TRUE, rotation = "none", nstart = 100, 
smartStart = NULL, seed = 1234)
}

\arguments{
  \item{data}{Continuous dataset}
  \item{nclus}{Number of clusters}
  \item{ndim}{Dimensionality of the solution}
  \item{alpha}{Adjusts for the relative importance of RKM and FKM in the objective function; \code{alpha = 1} reduces to PCA, \code{alpha = 0.5} to reduced K-means, and \code{alpha = 0} to factorial K-means}
  \item{method}{Specifies the method. Options are RKM for reduced K-means and FKM for factorial K-means (default = \code{"RKM"})}
  \item{center}{A logical value indicating whether the variables should be shifted to be zero centered (default = \code{TRUE)}}
  \item{scale}{A logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place (default = \code{TRUE)}}
  \item{rotation}{Specifies the method used to rotate the factors. Options are \code{none} for no rotation, \code{varimax} for varimax rotation with Kaiser normalization and \code{promax} for promax rotation (default = \code{"none"})}
  \item{nstart}{Number of starts (default = 100)}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when smartStart = NULL. The default value is 1234}
}


\value{
\item{obscoord}{Object scores}
\item{attcoord}{Variable scores}
\item{centroid}{Cluster centroids}
\item{cluID}{Cluster membership}
\item{criterion}{Optimal value of the objective function}
\item{size}{The number of objects in each cluster}
\item{scale}{A copy of \code{scale} in the return object}
\item{center}{A copy of \code{center} in the return object}
\item{nstart}{A copy of \code{nstart} in the return object}
\item{odata}{A copy of \code{data} in the return object}
}
\references{
De Soete, G., and Carroll, J. D. (1994). K-means clustering in a low-dimensional Euclidean space. In Diday E. et al. (Eds.), \emph{New Approaches in Classification and Data Analysis}, Heidelberg: Springer, 212-219.

Vichi, M., and Kiers, H.A.L. (2001). Factorial K-means analysis for two-way data. \emph{Computational Statistics and Data Analysis}, 37, 49-64.
}

\seealso{
\code{\link{clusmca}}, \code{\link{clusval}}, \code{\link{tuneclus}}
}
\examples{
#Reduced K-means with 3 clusters in 2 dimensions after 10 random starts
data(macro)
outRKM = cluspca(macro, 3, 2, method = "RKM", rotation = "varimax", scale = FALSE, nstart = 10)
outRKM
#Scatterplot (dimensions 1 and 2) and cluster description plot
plot(outRKM, cludesc = TRUE)

#Factorial K-means with 3 clusters in 2 dimensions 
#with a Reduced K-means starting solution
data(macro)
outFKM = cluspca(macro, 3, 2, method = "FKM", rotation = "varimax", 
scale = FALSE, smartStart=outRKM$cluID)
outFKM
#Scatterplot (dimensions 1 and 2) and cluster description plot
plot(outFKM, cludesc = TRUE)

#To get the Tandem approach (PCA + K-means)
outTandem = cluspca(macro, 3, 2, alpha = 1, scale = FALSE)
plot(outTandem)

#Cluster quality assessment in the low-dimensional space
clusval(outRKM, dst = "low")
clusval(outFKM, dst = "low")
clusval(outTandem, dst = "low")
}
