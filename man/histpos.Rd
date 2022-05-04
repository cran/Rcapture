\name{histpos}
\alias{histpos}
\alias{histpos.t}
\alias{histpos.0}
\title{ Possible Capture Histories  }
\description{ 
\code{histpos.t} builds the matrix of observable capture histories in terms of captures and misses for a capture-recapture 
experiment with \code{t} capture occasions.

\code{histpos.0} builds the matrix of observable capture histories in terms of number of captures for each primary period of 
a robust design with \code{vt} capture occasions. For closed populations, \code{vt=t} and \code{histpos.0} simply returns \code{t:1}.
}
\usage{
histpos.t(t)
histpos.0(vt)
}
\arguments{
  \item{t}{ The number of capture occasions. }
  \item{vt}{ A vector containing the numbers of capture occasions for each primary sampling period of a robust design. 
             The length of this vector equals the number of primary sampling periods (noted \eqn{I}). }
}
\value{
\code{histpos.t} gives a \eqn{2^t-1} by \eqn{t} matrix with rows representing capture histories. This matrix contains only 
zeros and ones.

\code{histpos.0} gives a \eqn{(t_1+1)\times(t_2+1)\times \ldots \times(t_I+1)-1}{(t1+1)*(t2+1)* \ldots *(tI+1)-1} by 
\eqn{I} matrix with rows representing possible capture histories in terms of number of captures. This matrix elements are 
integers between 0 and \eqn{t_i}{ti} inclusively. 
}
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \doi{10.18637/jss.v019.i05}.
}
\author{ Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and Sophie Baillargeon }

\note{ 
The \code{histpos.t} function is called by \code{\link{descriptive}}, \code{\link{closedp}}, \code{\link{closedp.bc}}, \code{\link{closedp.Mtb}}, \code{\link{openp}}, \code{\link{robustd.t}}.

The \code{histpos.0} function is called by \code{\link{robustd.0}}.
}
\seealso{ \code{\link{periodhist}}, \code{\link{descriptive}} }
\examples{
histpos.t(5)
histpos.0(5)

histpos.0(rep(5, 3))
}
\keyword{manip}
