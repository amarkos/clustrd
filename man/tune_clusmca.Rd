\name{tune_clusmca}
\alias{tune_clusmca}
\title{
Methods for categorical data with cluster quality assessment.
}
\description{
This function facilitates the selection of the appropriate number of clusters and dimensions for joint dimension reduction and clustering of categorical data. 
}
\usage{
tune_clusmca(data, nclusrange = 2:7, ndimrange = 2:4, method = "clusCA", 
criterion = "asw", dst = "full", alpha = .5, nstart = 100, 
smartStart = NULL, seed = 1234)
}
\arguments{
  \item{data}{Categorical dataset}
  \item{nclusrange}{An integer vector with the range of numbers of clusters which are to be compared by the cluster validity criteria. Note: the number of clusters should be greater than one}
  \item{ndimrange}{An integer vector with the range of dimensions which are to be compared by a cluster quality criterion}
  \item{criterion}{One of \code{asw}, \code{ch} or \code{crit}. Determines whether average silhouette width, Calinski-Harabasz index or objective value of the selected method is used \code{(default = "asw")}}
  \item{dst}{Specifies the data used to compute the distances between objects. Options are \code{full} for the original data (after possible scaling) and \code{low} for the object scores in the low-dimensional space \code{(default = "full")}}
  \item{method}{Specifies the method. Options are MCAk for MCA K-means, iFCB for Iterative Factorial Clustering of Binary variables and clusCA for Cluster Correspondence Analysis (default = \code{clusCA)}.}
  \item{alpha}{Non-negative scalar to adjust for the relative importance of of MCA and K-means in the solution \code{(default = .5)}. Works only in combination with \code{method = "MCAk"}.}
  \item{nstart}{Number of random starts.}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when \code{smartStart = NULL}. The default value is 1234.
}
}

\value{
\item{clusmcaobj}{The output of the optimal run of the \code{clusmca()} function}
\item{nclusbest}{The optimal number of cluster}
\item{ndimbest}{The optimal number of dimensions}
\item{critbest}{The optimal criterion value for \code{nclusbest} clusters and \code{ndimbest} dimensions}
\item{critgrid}{Matrix of size \code{nclusrange x ndimrange} with criterion values for the specified range of numbers of clusters and numbers of dimensions (values are calculated for the number of clusters greater or equal to the number of dimensions; otherwise values are left blank)}
}

\seealso{
\code{\link{clusmca}}, \code{\link{tune_cluspca}}
}

\examples{
data(underwear)
bestclusCA = tune_clusmca(underwear[,2:3], 3:4, 2, method = "clusCA", criterion = "asw")
# plot(bestclusCA$clusmcaobj)
}

