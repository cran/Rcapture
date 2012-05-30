\name{rbvole}
\alias{rbvole}
\docType{data}
\title{ Robust Design Data for Red-Back Voles }
\description{
  This data set contains robust design capture history data for red-back voles. 
}
\usage{data(rbvole)}
\format{
66 by 19 numeric matrix, with the following columns:
  \describe{
    \item{\code{c11}, \code{c12}, \code{c13}}{ Capture histories for the three capture occasions within primary period 1}
    \item{\code{c21}, \code{c22}, \code{c23}}{ Capture histories for the three capture occasions within primary period 2}
    \item{\code{c31}, \code{c32}, \code{c33}}{ Capture histories for the three capture occasions within primary period 3}
    \item{\code{c41}, \code{c42}, \code{c43}}{ Capture histories for the three capture occasions within primary period 4}
    \item{\code{c51}, \code{c52}, \code{c53}}{ Capture histories for the three capture occasions within primary period 5}
    \item{\code{c61}, \code{c62}, \code{c63}}{ Capture histories for the three capture occasions within primary period 6}
    \item{freq}{ Observed frequencies for each capture history}
    }
}
\details{
  Data collection was carried out by Etcheverry and al.. The capture occasions represent three consecutive days of trapping in 
  May 1999, July 1999, August 1999, May 2000, July 2000 and August 2000 in the Duchenier conservation area in southeastern Quebec. 
  This data set is analyzed in Rivest and Daigle (2004).

  Each row of this data set represents an observed capture history followed by its frequency. 
}
\references{
   Rivest, L.P. and Daigle, G. (2004) Loglinear models for the robust design in mark-recapture experiments. \emph{Biometrics}, \bold{60}, 100--107.
}
\examples{
  # According to Rivest and Daigle (2004), a good robust design model
  # for this data set is formed of an Mth Chao model for each period.
  # This model can be fitted as follows.
  
### data(rbvole)
### memory.limit(size=2000)
### rd<-robustd.t(rbvole,dfreq=TRUE,vt=rep(3,6),vm="Mth",vh="Chao")

  # WARNING : Because the data has 18 capture occasions, fitting this
  # model uses a lot of memory. Its runtime is several minutes long.
}
\keyword{datasets}
