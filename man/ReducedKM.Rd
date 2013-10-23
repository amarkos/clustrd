\name{ReducedKM}
\alias{ReducedKM}
\title{ReducedKM}
\description{Implements Reduced k-means (De Soete and Carroll, 1994) which combines k-means for clustering with PCA for dimension reduction.}
\usage{ReducedKM(data, nclus, ndim, nstart = 100, smartStart = F)}
\arguments{
  \item{data}{quantitative dataset}
  \item{nclus}{number of clusters}
  \item{ndim}{dimensionality of the solution}
  \item{nstart}{number of starts}
  \item{smartStart}{If TRUE then starting values are obtained with k-means}
}
\value{
  \item{obscoord}{object scores}
  \item{attcoord}{variable loadings}
  \item{centroid}{cluster centroids}
  \item{cluID}{cluster membership}
  \item{criterion}{optimal value of the objective function}
}
\references{De Soete, G. and Carroll, J. D. (1994). K-means clustering in a low-dimensional Euclidean space. In Diday E. et al. (Eds.), New Approaches in Classification and Data Analysis, Heidelberg: Springer, 212-219.}
\author{Markos, A. \email{amarkos@gmail.com}, Iodice D'Enza, A. \email{iodicede@gmail.com} and Van de Velden, M. \email{vandevelden@ese.eur.nl}}
\seealso{\code{\link{FactorialKM}}}
\examples{
   data(macro)
   outr <- ReducedKM(macro,3,2,nstart=1,smartStart=TRUE)
   plotrd(outr,what=c("all","none"),obslabel=rownames(macro),density=FALSE)
}
\keyword{cluster}
\keyword{clustering}
