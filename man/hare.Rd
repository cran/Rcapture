\name{hare}
\alias{hare}
\docType{data}
\title{  Snowshoe Hare Data }
\description{
  This data set contains capture-recapture data for snowshoe hares. 
}
\usage{hare}
\format{
68 by 6 numeric matrix, with the following columns:
  \describe{
    \item{\code{c1}, \code{c2}, \code{c3}, \code{c4}, \code{c5}, \code{c6}}{Capture histories from the six capture occasions}
  }
}
\details{
  This data set is analyzed in Cormack (1989) and Agresti (1994).

  Each row of this data set represents the capture history of one animal. 
}
\references{
Agresti, A.  (1994) Simple capture-recapture models permitting unequal catchability and variable sampling effort. \emph{Biometrics}, \bold{50}, 494--500.

Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \doi{10.18637/jss.v019.i05}.

Cormack, R. M. (1989) Loglinear models for capture-recapture. \emph{Biometrics}, \bold{45}, 395--413.
}
\examples{
desc <- descriptive(hare)
plot(desc)
  # The fi plot shows that the two animals caught on all occasions create 
  # some heterogeneity in the capture probabilities.

closedp(hare)
  # The best fitting model Mth Poisson2(N = 81.1, s.e. = 5.7) has an AIC of 146. 
closedpCI.t(hare, m = "Mth", h = "Poisson", h.control = list(theta = 2))

  # One can compare the fit of this model with that obtained by removing the
  # 2 hares caught 6 times. This can be done by adding a column to the design 
  # matrix for Mt taking the value 1 for the capture history (1,1,1,1,1,1).
col <- rep(0, 2^6-1)
mat <- histpos.t(6)
col[rowSums(mat) == 6] <- 1
closedpCI.t(hare, mX = cbind(mat, col), mname = "Mt without 111111")
  # This gives N = 76.8 (s.e. = 3.9) with an AIC of 146.
}
\keyword{datasets}
