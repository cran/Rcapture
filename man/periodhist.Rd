\name{periodhist}
\alias{periodhist}
\title{  Merge capture occasions }
\description{
This function produces a reduced matrix of capture histories from a complete one by merging together some capture occasions.

It can also be used to change the format of a capture-recapture data set with complete capture histories: it transforms a data 
set with one row per captured unit to a data set with one row per capture history followed by its frequency.
}
\usage{
periodhist(X, dfreq=FALSE, vt, drop=TRUE)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} 
            for a description of the accepted formats).  }
  \item{dfreq}{ A logical. By default FALSE, which means that \code{X} has one row per unit. If TRUE, it 
                indicates that the matrix \code{X} contains frequencies in its last column.}
  \item{vt}{ A vector containing the numbers of capture occasions for each pooled capture occasions. The length of 
             this vector equals the number of capture occasions in the reduced matrix (noted \eqn{I}). }
  \item{drop}{ A logical, by default TRUE, meaning that in the output matrix the unobserved 
               capture histories (frequency of 0) are omitted. To keep them, 
               \code{drop} must be set to FALSE.}
}
\details{
  This function is useful when using an open population model to analyze a robust design data set. It can be used to 
  reduce the data set to one observation per primary period. A one represents a unit caught at least once during the period 
  and a zero a unit never caught.
  
  It is also useful for experiments with a large number of capture occasions but a limited number of catches,  
  especially when there is no capture on some occasions. In such cases, one can pool together some capture occasions 
  with the \code{periodhist} function.  
  
  If a data matrix produced by the \code{periodhist} function is given as an argument to a \pkg{Rcapture} function, 
  the \code{dfreq} argument must be set to \code{TRUE}.
}
\value{
 A \eqn{2^I-1} by \eqn{I+1} matrix of all the pooled capture histories, 
 with their observed frequencies. If \code{drop=TRUE} (the default), the number of rows
 of this matrix is in fact \eqn{2^I-1} minus the number of unobserved capture histories.
 }
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.
}
\author{ Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and Sophie Baillargeon }

\seealso{ \code{\link{openp}} }
\examples{
data(mvole)
mvole.op <- periodhist(mvole, vt = rep(5, 6))
openp(mvole.op, dfreq = TRUE)
}
\keyword{ manip }
