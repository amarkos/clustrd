\name{underwear}
\alias{underwear}
\docType{data}
\title{South Korean Underwear}
\description{
The dataset comes from a large survey conducted by a South Korean underwear manufacturer in 1997. 664 South Korean consumers were asked to provide responses for three multiple-choice items: attributes when considering a brand of underwear to purchase (15 attributes), preferred brand of underwear (8 brands) and consumer age (3 levels).
}
\usage{data(underwear)}
\format{
  A data frame with 664 observations on the following variables.
  \describe{
    \item{\code{brand}}{categorical: \code{1. BYC}, \code{2. TRY}, \code{3. VICMAN}, \code{4. James Dean}, \code{5. Michiko-London}, \code{6. Benetton}, \code{7. Bodyguard}, \code{8. Calvin Klein}}
    \item{\code{atts}}{categorical: \code{1. Comfortable}, \code{2. Smooth}, \code{3. Superior fabrics}, \code{4. Reasonable price}, \code{5. Fashionable design}, \code{6. Favorable advertisements}, \code{7. Trendy color}, \code{8. Good design}, \code{9. Various colors}, \code{10. Elastic}, \code{11. Store is near}, \code{12. Excellent fit}, \code{13. Design quality}, \code{14. Youth appeal},
\code{15. Various sizes}}
    \item{\code{age}}{categorical: \code{1. 10-29}, \code{2. 30-49}, \code{3. 50 and over}}
  }
}

\source{
Hwang, H., Dillon, W. R. & Takane, Y. (2006). An extension of multiple correspondence analysis for identifying heterogenous subgroups of respondents. \emph{Psychometrika}, 71, 161-171.
}
\examples{
data(underwear)
}
\keyword{datasets}


