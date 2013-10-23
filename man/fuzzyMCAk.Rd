\name{fuzzyMCAk}
\alias{fuzzyMCAk}
\title{fuzzyMCAk}
\description{Implements Fuzzy Cluster MCA (Hwang, Dillon and Takane, 2010) which combines Multiple Correspondence Analysis for dimension reduction with fuzzy c-means (Bezdek, 1981) for clustering.}
\usage{
fuzzyMCAk(data,nclus=3,ndim=2,nstart=1)}
\arguments{
\item{data}{categorical dataset}                                                    
  \item{nclus}{number of clusters}
  \item{ndim}{dimensionality of the solution}                           
  \item{nstart}{number of random starts}                                
  }
\value{
  \item{obscoord}{object scores}
  \item{attcoord}{attribute scores}
  \item{centroid}{cluster centroids}
  \item{cluID}{hard cluster membership}
  \item{U}{fuzzy cluster membership} 
  \item{FPI}{Fuzziness Performance Index}
  \item{MPE}{Modified Partition Entropy}
  \item{criterion}{optimal value of the objective function}
}
\references{Bezdek, J. C. (1981). Pattern recognition with fuzzy objective function algorithms. Kluwer Academic Publishers.

Hwang, H., Dillon, W. R. and Takane, Y. (2010). Fuzzy Cluster Multiple Correspondence Analysis. Behaviormetrika, 37(2), 111-133.
}
\author{Markos, A. \email{amarkos@gmail.com}, Iodice D'Enza , A. \email{iodicede@gmail.com} and Van de Velden, M. \email{vandevelden@ese.eur.nl}}
\seealso{\code{\link{groupals}}, \code{\link{iFCB}}, \code{\link{MCAk}}}
\examples{
     data(underwear)  
     outfMCAk <- fuzzyMCAk(underwear[c(1:200),c(1:2)],nclus=3,ndim=2,nstart=1)
     plotrd(outfMCAk)
     }
\keyword{cluster}
\keyword{clustering}% __ONLY ONE__ keyword per line
