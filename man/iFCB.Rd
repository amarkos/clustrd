\name{iFCB}
\alias{iFCB}
\title{iFCB}
\description{Implements iFCB (Iodice D'Enza and Palumbo, 2013) which combines Nonsymmetric Correspondence Analysis for dimension reduction with k-means for clustering.}
\usage{iFCB(data,nclus,ndim,nstart=100,smartStart=F,seed=1234)}
\arguments{
  \item{data}{categorical dataset}                                                    
  \item{nclus}{number of clusters}
  \item{ndim}{dimensionality of the solution}                           
  \item{nstart}{number of random starts}                                
  \item{smartStart}{If TRUE then starting values are obtained with k-means}
  \item{seed}{seed is used to set the random number seed when smartStart = FALSE}
  }
\value{
  \item{obscoord}{object scores}
  \item{attcoord}{attribute scores}
  \item{centroid}{cluster centroids}
  \item{cluID}{cluster membership}
  \item{criterion}{optimal value of the objective function}
}
\references{Iodice D' Enza, A. and Palumbo, F. (2013). Iterative factor clustering of binary data. Computational Statistics, 28 (2), 789-807.}
\author{Markos, A. \email{amarkos@gmail.com}, Iodice D'Enza, A. \email{iodicede@gmail.com} and Van de Velden, M. \email{vandevelden@ese.eur.nl}}
\seealso{\code{\link{MCAk}}, \code{\link{fuzzyMCAk}}, \code{\link{groupals}}}
\examples{
     data(underwear)
     attlab = c(c(1:15),"by","tr","vm","jd","ml","bn","bg","ck","a1","a2","a3")
     outiFCB <- iFCB(underwear,nclus=3,ndim=2,nstart=1,smartStart=TRUE,seed=1234)
     plotrd(outiFCB,attlabel=attlab)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{clustering}% __ONLY ONE__ keyword per line
