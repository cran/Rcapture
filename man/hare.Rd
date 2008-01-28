\name{hare}
\alias{hare}
\docType{data}
\title{  Closed Population Data for Snowshoe Hares }
\description{
  This data set contains closed population capture history data for snowshoe hares. 
}
\usage{data(hare)}
\format{
  \describe{
    \item{\code{c1}, \code{c2}, \code{c3}, \code{c4}, \code{c5}, \code{c6}}{Capture histories for the six capture occasions}
  }
}
\details{
  This data set is analysed in Cormack (1989) and Agresti (1994).

  This data set's format is the default one, i.e. each row represents the capture history of one animal. 
}
\references{
Agresti, A.  (1994). Simple capture-recapture models permitting unequal catchability and variable sampling effort. \emph{Biometrics}, 50, 494--500.

Baillargeon, S. and Rivest, L.P. (2007). Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), http://www.jstatsoft.org/

Cormack, R. M. (1989). Log-linear models for capture-recapture. \emph{Biometrics}, \bold{45}, 395--413.
}
\examples{
data(hare)
desc<-descriptive(hare)
plot(desc)
  # The fi plot shows that the two animals caught on all occasions create 
  # some heterogeneity in the capture probabilities.

closedp(hare)
  # The best fitting model Mth Poisson2(N = 81.1, s.e.=5.7) has an AIC of 146. 

  # One can compare the fit of this model with that obtained by removing the
  # 2 hares caught 6 times. This can be done by adding a column to the design 
  # matrix for Mt taking the value 1 for the capture history (1,1,1,1,1,1).
col<-rep(0,2^6-1)
mat<-histpos.t(6)
col[apply(mat,1,sum)==6]<-1
closedp.mX(hare,mX=cbind(mat,col),mname="Mt without 111111")
  # This gives N = 76.8 (s.e.=3.9) with an AIC of 146.

  # Profile likelihood confidence intervals for both models:
profileCI(hare,m="Mth",h="Poisson",a=2)
profileCI(hare,mX=cbind(mat,col),mname="Mt without 111111")
}
\keyword{datasets}
