\name{lesbian}
\alias{lesbian}
\docType{data}
\title{ Epidemiological capture-recapture Data on the lesbian population}
\description{
  Epidemiological capture-recapture data on the lesbian population from four four organisations
that serve the lesbian and gay population of Allegheny County, Pennsylvania, United States and maintain large mailing lists.
}
\usage{lesbian}
\format{
15 by 5 numeric matrix, with the following columns:
  \describe{
    \item{\code{A}, \code{B}, \code{C}, \code{D}}{ Capture histories from the four organisations mailling lists}
    \item{\code{freq}}{ Observed frequencies for each capture history}
    }
}
\details{
  The capture histories were obtained by matching the names on the four mailling lists. 

  Each row of this data set represents an observed capture history followed by its frequency. 
}
\source{
Aaron, D.J., Chang, Y.-F., Markovic, N., LaPorte, R. E. (2003) Estimating the lesbian population: a capture-recapture approach. \emph{Journal of Epidemiology and Community Health }, \bold{57}(3), 207--209.
}
\examples{
closedp(lesbian, dfreq = TRUE)

# According to the BIC, the best model is Mth Darroch.
# Let's see if adding  interactions between capture 
# histories to the model could improve the model's fit.

closedpMS.t(lesbian, dfreq = TRUE, h = "Darroch")

# According to the BIC, the best heterogeneous Darroch model
# contains the double interactions 12, 13, 14. 
}
\keyword{datasets}
