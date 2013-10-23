\name{plotrd}
\alias{plotrd}
\title{plotrd}
\description{Plot function that creates a ggplot2 based density map of the object scores and a scatter plot of both the attribute scores and the centroids. Creates a .png file in the working directory.}
\usage{plotrd(clustrdOut,what=c("all","all"),obslabel=0,attlabel=0,density=T,fname=0)}
\arguments{
  \item{clustrdOut}{list object created by one of the methods.}
  \item{what}{Vector of two character strings specifying the contents of the plot. First entry sets the rows and the second entry the columns. Allowed values are "all" (all available points, default) "none" (no points are displayed) }
  \item{obslabel}{the object labels; if not provided, no labeling is applied.}
  \item{attlabel}{the attribute labels; if not provided, default labeling is applied.}
  \item{density}{if TRUE a density map of the objects is created.}
  \item{fname}{name of the png file to be created. Default file name is "figure".}
}
\value{a ggplot object}
\author{Markos, A. \email{amarkos@gmail.com}, Iodice D'Enza, A. \email{iodicede@gmail.com} and Van de Velden, M. \email{vandevelden@ese.eur.nl}}
\keyword{cluster}
\keyword{clustering}% __ONLY ONE__ keyword per line


