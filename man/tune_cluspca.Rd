\name{tune_cluspca}
\alias{tune_cluspca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Methods for continuous data with cluster quality assessment.
}
\description{
This function facilitates the selection of the appropriate number of clusters and dimensions for joint dimension reduction and clustering of continuous data. 
}
\usage{
tune_cluspca(data, nclusrange = 2:7, ndimrange = 2:4, criterion = "asw", 
dst = "full", alpha = NULL, method = "RKM", center = TRUE, scale = TRUE, 
rotation = "none", nstart = 10, smartStart = NULL, seed = 1234)
}

\arguments{
  \item{data}{Continuous dataset}
  \item{nclusrange}{An integer vector with the range of numbers of clusters which are to be compared by the cluster validity criteria}
  \item{ndimrange}{An integer vector with the range of dimensions which are to be compared by the cluster validity criteria}
  \item{criterion}{One of \code{asw}, \code{ch} or \code{crit}. Determines whether average silhouette width, Calinski-Harabasz index or objective value of the selected method is used \code{(default = "asw")}}
  \item{dst}{Specifies the data used to compute the distances between objects. Options are \code{full} for the original data (after possible scaling) and \code{low} for the object scores in the low-dimensional space \code{(default = "full")}}
    \item{alpha}{Adjusts for the relative importance of the two terms of Clustering and Dimension Reduction; \code{alpha = 1} reduces to PCA, \code{alpha = 0.5} to reduced K-means, and \code{alpha = 0} to factorial K-means}
  \item{method}{Specifies the method. Options are RKM for reduced K-means and FKM for factorial K-means (\code{default = "RKM"}).}
  \item{center}{A logical value indicating whether the variables should be shifted to be zero centered \code{(default = TRUE)}}
  \item{scale}{A logical value indicating whether the variables should be scaled to have unit variance before the analysis takes place \code{(default = TRUE)}}
  \item{rotation}{Specifies the method used to rotate the factors. Options are none for no rotation, varimax for varimax rotaion with Kaiser normalization and promax
for promax rotation \code{(default = "none")}}
  \item{nstart}{Number of starts}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when smartStart = NULL. The default value is 1234}
}

\value{
\item{cluspcaobj}{The output of the optimal run of the \code{cluspca()} function}
\item{nclusbest}{The optimal number of cluster}
\item{ndimbest}{The optimal number of dimensions}
\item{critbest}{The optimal criterion value for \code{nclusbest} clusters and \code{ndimbest} dimensions}
\item{critgrid}{Matrix of size \code{nclusrange x ndimrange} with criterion values for the specified ranges of numbers of clusters and numbers of dimensions (values are calculated for the number of clusters greater than the number of dimensions; otherwise values are left blank)}
}


\seealso{
\code{\link{cluspca}}, \code{\link{tune_clusmca}}
}
\examples{
data(macro)
bestRKM = tune_cluspca(macro, 3:4, 2:3, method = "RKM", criterion = "asw", dst = "low")
plot(bestRKM$cluspcaobj)
}

