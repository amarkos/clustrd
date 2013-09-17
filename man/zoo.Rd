\name{zoo}
\alias{zoo}
\docType{data}
\title{zoo dataset}
\description{A dataset containing 16 binary and one ordinal variable on 101 animals.}
\usage{data(zoo)}
\format{
  A data frame with 101 animals on the following 17 attributes.
  \describe{
    \item{\code{hair}}{binary}
    \item{\code{feathers}}{binary}
    \item{\code{eggs}}{binary}
    \item{\code{milk}}{binary}
    \item{\code{airborne}}{binary}
    \item{\code{aquatic}}{binary}
    \item{\code{predator}}{binary}
    \item{\code{toothed}}{binary}
    \item{\code{backbone}}{binary}
    \item{\code{breathes}}{binary}
    \item{\code{venomous}}{binary}
    \item{\code{fins}}{binary}
    \item{\code{legs}}{numeric}
    \item{\code{tail}}{a numeric vector}
    \item{\code{domestic}}{binary}
    \item{\code{catsize}}{binary}
    \item{\code{truecluID}}{a factor with levels \code{mammal} \code{bird} \code{reptile} \code{fish} \code{amphibian} \code{insect} \code{mollusc.et.al}}
  }
}
\source{Newman, D.J. & Hettich, S. & Blake, C.L. & Merz, C.J. (1998). UCI Repository of machine learning databases [http://www.ics.uci.edu/~mlearn/MLRepository.html]. Irvine, CA: University of California, Department of Information and Computer Science.}
\examples{
data(zoo)
summary(zoo)
}
\keyword{datasets}
