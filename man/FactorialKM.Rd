\name{FactorialKM}
\alias{FactorialKM}
\title{FactorialKM}
\description{Implements Factorial k-means (Vichi and Kiers, 2001) which combines Principal Component Analysis for dimension reduction with k-means for clustering.}
\usage{FactorialKM(data, nclus, ndim, nstart = 100, smartStart = F)}
\arguments{
  \item{data}{quantitative dataset}
  \item{nclus}{number of clusters}
  \item{ndim}{dimensionality of the solution}
  \item{nstart}{number of starts}
  \item{smartStart}{If TRUE then starting values are obtained with k-means}
}
\value{
  \item{obscoord}{object scores}
  \item{attcoord}{attribute scores}
  \item{centroid}{cluster centroids}
  \item{cluID}{cluster membership}
  \item{criterion}{optimal value of the objective function}
}
\references{Vichi, M. and Kiers, H.A.L. (2001). Factorial k-means analysis for two-way data. Computational Statistics and Data Analysis, 37, 49-64.}
\author{Markos, A. \email{amarkos@gmail.com}, Iodice D'Enza , A. \email{iodicede@gmail.com} and Van de Velden, M. \email{vandevelden@ese.eur.nl}}
\seealso{\code{\link{ReducedKM}}}
\examples{
   data(macro)
   outf <- FactorialKM(macro,3,2,nstart=1,smartStart=TRUE)
   plotrd(outf,what=c("all","none"),obslabel=rownames(macro),density=FALSE)
}
\keyword{cluster}
\keyword{clustering}% __ONLY ONE__ keyword per line
