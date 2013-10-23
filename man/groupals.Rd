\name{groupals}
\alias{groupals}
\title{groupals}
\description{Implements a constrained homogeneity analysis to take into account the cluster structure of objects as described in Van Buuren and Heiser (1989).}
\usage{groupals(data,nclus,ndim,nstart=100,smartStart=F,seed=1234)}
\arguments{
\item{data}{categorical dataset}                                                    
  \item{nclus}{number of clusters}
  \item{ndim}{dimensionality of the solution}                           
  \item{nstart}{number of random starts}                                
  \item{smartStart}{If TRUE then starting values are obtained with fuzzy c-means}
  \item{seed}{seed is used to set the random number seed when smartStart = FALSE}
  }
\value{
  \item{obscoord}{object scores}
  \item{attcoord}{attribute scores}
  \item{centroid}{cluster centroids}
  \item{cluID}{cluster membership}
  \item{criterion}{optimal value of the objective function}
}
\references{Van Buuren, S. and Heiser, W. J. (1989). Clustering n objects in k groups under optimal scaling of variables, Psychometrika, 54, 699-706.}
\author{Markos, A. \email{amarkos@gmail.com}, Iodice D'Enza, A. \email{iodicede@gmail.com} and Van de Velden, M. \email{vandevelden@ese.eur.nl}}
\seealso{\code{\link{MCAk}}, \code{\link{fuzzyMCAk}}, \code{\link{iFCB}}}
\examples{
     data(underwear)
     attlab = c(c(1:15),"by","tr","vm","jd","ml","bn","bg","ck","a1","a2","a3")
     outgroupals <- groupals(underwear,nclus=3,ndim=2,nstart=1,smartStart=TRUE,seed=1234)
     plotrd(outgroupals,attlabel=attlab)
}

\keyword{cluster}
\keyword{clustering}

