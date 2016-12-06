\name{clusmca}
\alias{clusmca}
\title{
Joint  dimension  reduction  and  clustering  of  categorical  data.
}

\description{
This function implements MCA K-means (Hwang, Dillon and Takane, 2006), i-FCB (Iodice D' Enza and Palumbo, 2013) and Cluster Correspondence Analysis (van de Velden, Iodice D' Enza and Palumbo, 2016). The methods combine variants of Correspondence Analysis for dimension reduction with K-means for clustering.
}
\usage{
clusmca(data, nclus, ndim, method = "clusCA", alpha = .5,
nstart = 10, smartStart = NULL, gamma = TRUE, seed = 1234)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Categorical dataset}
  \item{nclus}{Number of clusters}
  \item{ndim}{Dimensionality of the solution}
  \item{method}{Specifies the method. Options are MCAk for MCA K-means, iFCB for Iterative Factorial Clustering of Binary variables and clusCA for Cluster Correspondence Analysis (default = \code{"clusCA")}}
  \item{alpha}{Non-negative scalar to adjust for the relative importance of of MCA and K-means in the solution (default = .5). Works only in combination with \code{method = "MCAk"}}
  \item{nstart}{Number of random starts}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{gamma}{Scaling parameter that leads to a similar spread in the object and attribute points (default = \code{TRUE})}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when \code{smartStart = NULL}. The default value is 1234}
}

\value{
\item{obscoord}{Object scores}
\item{attcoord}{Varable scores}
\item{centroid}{Cluster centroids}
\item{cluID}{Cluster membership}
\item{criterion}{Optimal value of the objective criterion}
\item{csize}{Cluster size}
\item{nstart}{A copy of \code{nstart} in the return object}
\item{odata}{A copy of \code{data} in the return object}
}
\references{
Hwang, H., Dillon, W. R. and Takane, Y. (2006). An extension of multiple correspondence analysis for identifying heterogenous subgroups of respondents, \emph{Psychometrika}, 71, 161-171.

Iodice D'Enza, A. and Palumbo, F. (2013). Iterative factor clustering of binary data. \emph{Computational Statistics}, \emph{28}(2), 789-807.

van de Velden M., Iodice D'Enza, A. and Palumbo, F. (2016).  Cluster correspondence analysis.\emph{Psychometrika} (in press) DOI: 10.1007/s11336-016-9514-0
}

\seealso{
\code{\link{cluspca}}, \code{\link{tune_clusmca}}
}
\examples{
data(cmc)
# values of wife's age and number of children were categorized 
# into three groups based on quartiles
cmc$W_AGE = ordered(cut(cmc$W_AGE, c(16,26,39,49), include.lowest = TRUE))
levels(cmc$W_AGE) = c("16-26","27-39","40-49") 
cmc$NCHILD = ordered(cut(cmc$NCHILD, c(0,1,4,17), right = FALSE))
levels(cmc$NCHILD) = c("0","1-4","5 and above")
outclusMCA = clusmca(cmc[,-c(3,6,7)], 3, 2, method = "clusCA")
}
