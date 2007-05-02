\name{periodhist}
\alias{periodhist}
\title{  Capture Histories Between Primary Sampling Periods }
\description{
 This function produces the matrix of capture histories between primary sampling periods from a complete matrix of capture histories for a robust design.
}
\usage{
periodhist(X, dfreq=FALSE, vt)
}
\arguments{
  \item{X}{ The table of the observed capture histories in one of the two accepted formats. In the default format, it has one row per unit captured in the experiment and sum(vt) columns. In the alternative format, it contains one row per capture history followed by its frequency. In that case, \code{X} has sum(vt)+1 columns. The first sum(vt) columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. }
  \item{dfreq}{ This argument specifies the format of the data matrix \code{X}. By default, it is set to FALSE, which means that \code{X} has one row per unit. If it is set to TRUE, then the matrix \code{X} contains frequencies in its last column.}
  \item{vt}{ A vector containing the numbers of capture occasions for each primary sampling period of a robust design. The length of this vector equals the number of primary sampling periods (noted \eqn{I}). }
}
\details{
  This function is useful when using an open population model to analyse a robust design data set. It reduces the data set to one observation per primary period. The observation is one for a unit caught at least once during the period and zero for a unit never caught. If a data matrix produced by the \code{periodhist} function is given in argument to an \code{openp} function, the \code{dfreq} argument must be set to TRUE.
}
\value{
 A \eqn{2^I-1} by \eqn{I+1} matrix of all the possible between periods capture histories, with their observed frequencies. This data matrix has the alternative format.
 }
\references{
Baillargeon, S. and Rivest, L.P. (2007). The Rcapture package: Loglinear models for capture-recapture in R. \emph{Journal of Statistical Software}, to appear (available online at \url{http://www.jstatsoft.org/}).
}
\author{Sophie Baillargeon \email{sbaillar@mat.ulaval.ca} and Louis-Paul Rivest \email{lpr@mat.ulaval.ca} }
\seealso{ \code{\link{openp}} }
\examples{
data(mvole)
mvole.op<-periodhist(mvole,vt=rep(5,6))
openp(mvole.op, dfreq=TRUE)
}
\keyword{ manip }
