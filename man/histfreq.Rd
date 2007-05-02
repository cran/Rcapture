\name{histfreq}
\alias{histfreq.t}
\alias{histfreq.0}
\title{ Frequencies of the Possible Capture Histories }
\description{
 These functions compute the frequencies of all possible capture histories.\cr
 \code{histfreq.t} is used for the \eqn{2^t-1} observable capture histories genereted by \code{histpos.t} in terms of captures and misses for each of the \eqn{t} capture occasions. \cr
 \code{histfreq.0} is used for observable capture histories genereted by \code{histpos.0} in terms of the numbers of captures in the primary periods of a robust design with \code{vt} capture occasions.
}
\usage{
histfreq.t(X, dfreq=FALSE)
histfreq.0(X, dfreq=FALSE, vt)
}
\arguments{
  \item{X}{ The table of the observed capture histories in one of the two accepted formats. In the default format, it has one row per unit captured in the experiment. In this case, the number of columns in the table represents the number of capture occasions in the experiment (noted \eqn{t}). In the alternative format, it contains one row per capture history followed by its frequency. In that case, \code{X} has \eqn{t}+1 columns. The first \eqn{t} columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. }
  \item{dfreq}{ This argument specifies the format of the data matrix \code{X}. By default, it is set to FALSE, which means that \code{X} has one row per unit. If it is set to TRUE, then the matrix \code{X} contains frequencies in its last column.}
  \item{vt}{ A vector containing the numbers of capture occasions for each primary sampling period of a robust design. The length of this vector equals the number of primary sampling periods (noted \eqn{I}). }
}
\value{
\code{histfreq.t} returns a vector of frequencies for the \eqn{2^t - 1} or \eqn{2^{\mbox{sum(vt)}}-1}{2^sum(vt)-1} capture histories ordered as in \code{histpos.t}. 

\code{histfreq.0} returns a vector of frequencies for the \eqn{(t_1+1)\times(t_2+1)\times \cdots \times(t_I+1)-1}{(t1+1)*(t2+1)* \cdots *(tI+1)-1} capture histories ordered as in \code{histpos.0} .
}
\references{
Baillargeon, S. and Rivest, L.P. (2007). The Rcapture package: Loglinear models for capture-recapture in R. \emph{Journal of Statistical Software}, to appear (available online at \url{http://www.jstatsoft.org/}).
}
\author{Sophie Baillargeon \email{sbaillar@mat.ulaval.ca} and Louis-Paul Rivest \email{lpr@mat.ulaval.ca} }
\note{ 
The \code{histfreq.t} function is called by \code{\link{descriptive}}, \code{\link{closedp}}, \code{\link{closedp.bc}}, \code{\link{closedp.Mtb}}, \code{\link{closedp.mX}}, \code{\link{closedp.h}}, \code{\link{periodhist}}, \code{\link{openp}}, \code{\link{robustd.t}}, \code{\link{robustd.0}}.

The \code{histfreq.0} function is called by \code{\link{robustd.0}}. 
}
\seealso{ \code{\link{histpos.t}}, \code{\link{histpos.0}}, \code{\link{descriptive}} }
\examples{
data(mvole)

period3<-mvole[,11:15]
histfreq.t(period3)

histfreq.0(mvole[,1:15],vt=rep(5,3))
}
\keyword{manip}
