\name{tuneclus}
\alias{tuneclus}

\title{
Cluster quality assessment for a range of clusters and dimensions.
}
\description{
This function facilitates the selection of the appropriate number of clusters and dimensions for joint dimension reduction and clustering methods. 
}
\usage{
tuneclus(data, nclusrange = 2:7, ndimrange = 2:4, method = "RKM",  
criterion = "asw", dst = "full", alpha = NULL, alphak = NULL, 
center = TRUE, scale = TRUE, rotation = "none", nstart = 100, 
smartStart = NULL, seed = 1234)
}

\arguments{
  \item{data}{Continuous or Categorical dataset}
  \item{nclusrange}{An integer vector with the range of numbers of clusters which are to be compared by the cluster validity criteria. Note: the number of clusters should be greater than one}
  \item{ndimrange}{An integer vector with the range of dimensions which are to be compared by the cluster validity criteria}
  \item{method}{Specifies the method. Options are \code{RKM} for reduced K-means, \code{FKM} for factorial K-means, \code{MCAk} for MCA K-means, \code{iFCB} for Iterative Factorial Clustering of Binary variables and \code{clusCA} for Cluster Correspondence Analysis}
  \item{criterion}{One of \code{asw}, \code{ch} or \code{crit}. Determines whether average silhouette width, Calinski-Harabasz index or objective value of the selected method is used (default = \code{"asw")}}
  \item{dst}{Specifies the data used to compute the distances between objects. Options are \code{full} for the original data (after possible scaling) and \code{low} for the object scores in the low-dimensional space (default = \code{"full")}}
    \item{alpha}{Adjusts for the relative importance of RKM and FKM in the objective function; \code{alpha = 1} reduces to PCA, \code{alpha = 0.5} to reduced K-means, and \code{alpha = 0} to factorial K-means.}
    \item{alphak}{Non-negative scalar to adjust for the relative importance of MCA (\code{alphak = 1}) and K-means (\code{alphak = 0}) in the solution (default = .5). Works only in combination with \code{method = "MCAk"}}
  \item{center}{A logical value indicating whether the variables should be shifted to be zero centered (default = \code{TRUE)}}
  \item{scale}{A logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place (default = \code{TRUE)}}
  \item{rotation}{Specifies the method used to rotate the factors. Options are none for no rotation, varimax for varimax rotaion with Kaiser normalization and promax
for promax rotation (default = \code{"none")}}
  \item{nstart}{Number of starts (default = 100)}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when smartStart = NULL. The default value is 1234}
}

\value{
\item{clusobj}{The output of the optimal run of \code{cluspca()} or \code{clusmca()}}
\item{nclusbest}{The optimal number of clusters}
\item{ndimbest}{The optimal number of dimensions}
\item{critbest}{The optimal criterion value for \code{nclusbest} clusters and \code{ndimbest} dimensions}
\item{critgrid}{Matrix of size \code{nclusrange x ndimrange} with the criterion values for the specified ranges of clusters and dimensions (values are calculated only when the number of clusters is greater than the number of dimensions; otherwise values in the grid are left blank)}}

\references{Calinski, R.B., and Harabasz, J., (1974). A dendrite method for cluster analysis. \emph{Communications in Statistics}, 3, 1-27.

Kaufman, L., and Rousseeuw, P.J., (1990). \emph{Finding Groups in Data: An Introduction to Cluster Analysis}. Wiley, New York.}

\seealso{
\code{\link{cluspca}}, \code{\link{clusmca}}
}

\examples{
# Reduced K-means for a range of clusters and dimensions
# Cluster quality assessment based on the average silhouette width 
# in the low dimensional space
data(macro)
bestRKM = tuneclus(macro, 3:4, 2:3, method = "RKM", criterion = "asw", dst = "low",nstart = 10)
bestRKM
#plot(bestRKM$clusobj)

# Cluster Correspondence Analysis for a range of clusters and dimensions
# Cluster quality assessment based on the average silhouette width 
# in the full dimensional space
data(underwear)
bestclusCA = tuneclus(underwear, 3:4, 2:3, method = "clusCA", criterion = "asw", nstart = 10)
bestclusCA
#plot(bestclusCA$clusobj)
}
