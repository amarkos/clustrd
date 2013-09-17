\name{MCAk}
\alias{MCAk}
\title{MCAk}
\description{Implements Cluster MCA (Hwang, Dillon and Takane, 2006) which combines Multiple Correspondence Analysis for dimension reduction with k-means for clustering.}
\usage{
MCAk(data,nclus,ndim,nstart=100,smartStart=F)}
\arguments{

\item{data}{categorical dataset}                                                    
  \item{nclus}{number of clusters}
  \item{ndim}{dimensionality of the solution}                           
  \item{nstart}{number of random starts}                                
  \item{smartStart}{If TRUE then starting values are obtained with k-means}
  }
\value{
  \item{obscoord}{object scores}
  \item{attcoord}{attribute scores}
  \item{centroid}{cluster centroids}
  \item{cluID}{cluster membership}
  \item{criterion}{optimal value of the objective function}
}
\references{Hwang, H., Dillon, W. R. and Takane, Y. (2006). An extension of multiple correspondence analysis for identifying heterogenous subgroups of respondents, Psychometrika, 71, 161-171.}
\author{Markos, A. \email{<amarkos@gmail.com>}, Iodice D'Enza , A. \email{<iodicede@gmail.com>} and Van de Velden, M. \email{<vandevelden@ese.eur.nl>}}
\seealso{\code{\link{groupals}}, \code{\link{iFCB}}, \code{\link{fuzzyMCAk}}}
\examples{
     data(underwear)
     attlab = c(c(1:15),"by","tr","vm","jd","ml","bn","bg","ck","a1","a2","a3")
     outMCAk <- MCAk(underwear,nstart=1,nclus=3,ndim=2,smartStart=TRUE)
     plotrd(outMCAk,attlabel=attlab)
}
\keyword{cluster}
\keyword{clustering}% __ONLY ONE__ keyword per line
