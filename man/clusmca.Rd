\name{clusmca}
\alias{clusmca}
\title{
Joint  dimension  reduction  and  clustering  of  categorical  data.
}

\description{
This function implements MCA K-means (Hwang, Dillon and Takane, 2006), i-FCB (Iodice D' Enza and Palumbo, 2013) and Cluster Correspondence Analysis (van de Velden, Iodice D' Enza and Palumbo, 2017). The methods combine variants of Correspondence Analysis for dimension reduction with K-means for clustering.
}

\usage{
clusmca(data, nclus, ndim, method = "clusCA", alphak = .5,
nstart = 100, smartStart = NULL, gamma = TRUE, seed = 1234)
}

\arguments{
  \item{data}{Categorical dataset}
  \item{nclus}{Number of clusters}
  \item{ndim}{Dimensionality of the solution}
  \item{method}{Specifies the method. Options are MCAk for MCA K-means, iFCB for Iterative Factorial Clustering of Binary variables and clusCA for Cluster Correspondence Analysis (default = \code{"clusCA")}}
  \item{alphak}{Non-negative scalar to adjust for the relative importance of MCA (\code{alphak = 1}) and K-means (\code{alphak = 0}) in the solution (default = .5). Works only in combination with \code{method = "MCAk"}}
  \item{nstart}{Number of random starts (default = 100)}
  \item{smartStart}{If \code{NULL} then a random cluster membership vector is generated. Alternatively, a cluster membership vector can be provided as a starting solution}
  \item{gamma}{Scaling parameter that leads to similar spread in the object and variable scores (default = \code{TRUE})}
  \item{seed}{An integer that is used as argument by \code{set.seed()} for offsetting the random number generator when \code{smartStart = NULL}. The default value is 1234}
}

\value{
\item{obscoord}{Object scores}
\item{attcoord}{Variable scores}
\item{centroid}{Cluster centroids}
\item{cluID}{Cluster membership}
\item{criterion}{Optimal value of the objective criterion}
\item{size}{The number of objects in each cluster}
\item{nstart}{A copy of \code{nstart} in the return object}
\item{odata}{A copy of \code{data} in the return object}
}
\references{
Hwang, H., Dillon, W. R., and Takane, Y. (2006). An extension of multiple correspondence analysis for identifying heterogenous subgroups of respondents. \emph{Psychometrika}, 71, 161-171.

Iodice D'Enza, A., and Palumbo, F. (2013). Iterative factor clustering of binary data. \emph{Computational Statistics}, \emph{28}(2), 789-807.

van de Velden M., Iodice D' Enza, A., and Palumbo, F. (2017).  Cluster correspondence analysis. \emph{Psychometrika}, \emph{82}(1), 158-185.
}

\seealso{
\code{\link{cluspca}}, \code{\link{clusval}}, \code{\link{tuneclus}}
}
\examples{
data(cmc)
# Preprocessing: values of wife's age and number of children were categorized 
# into three groups based on quartiles
cmc$W_AGE = ordered(cut(cmc$W_AGE, c(16,26,39,49), include.lowest = TRUE))
levels(cmc$W_AGE) = c("16-26","27-39","40-49") 
cmc$NCHILD = ordered(cut(cmc$NCHILD, c(0,1,4,17), right = FALSE))
levels(cmc$NCHILD) = c("0","1-4","5 and above")

#Cluster Correspondence Analysis solution with 3 clusters in 2 dimensions after #10 random starts
outclusCA = clusmca(cmc, 3, 2, method = "clusCA", nstart = 10)
outclusCA
#Scatterplot (dimensions 1 and 2)
plot(outclusCA)

#MCA K-means solution with 3 clusters in 2 dimensions after 10 random starts
outMCAk = clusmca(cmc, 3, 2, method = "MCAk", nstart = 10)
outMCAk
#Scatterplot (dimensions 1 and 2)
plot(outMCAk)
}
