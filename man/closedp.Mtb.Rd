\name{closedp.Mtb}
\alias{closedp.Mtb}
\alias{print.closedp.Mtb}
\title{ Closed Population Capture-Recapture Model Mtb }
\description{
  This function fits model Mtb for closed populations in capture-recapture experiments.
}
\usage{
closedp.Mtb(X, dfreq=FALSE)

\method{print}{closedp.Mtb}(x, \dots)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} for a description of the accepted formats).  }
  \item{dfreq}{ A logical. By default FALSE, which means that \code{X} has one row per unit. If TRUE, it 
                indicates that the matrix \code{X} contains frequencies in its last column.}
  \item{x}{ An object, produced by the \code{closedp.Mtb} function, to print.}
  \item{\dots}{ Further arguments to be passed to \code{print.default}.}
  }
\details{ 
  The Mtb model is non-linear. It is fitted with the \code{optim} function instead of the \code{glm} function. 
  Therefore, the abundance estimate can be unstable. 

  For the model to be identifiable, the parameters are constrained in the following way: 
  \eqn{logit(c_i)=logit(p_i)+b}{logit(ci)=logit(pi)+b} for i in \eqn{2,\ldots,l}. 
}
\value{
\item{n}{ The number of captured units}
\item{results}{ A table containing the estimated population size and its standard error, the deviance, 
                the number of degrees of freedom and the Akaike's information criterion.}
\item{parMtb}{ Capture-recapture parameters estimates for model Mtb : the abundance N, \eqn{p_1}{p1} to \eqn{p_t}{pt},
               the probabilities of first capture for each capture occasion, and \eqn{c_2}{c2} to \eqn{c_t}{ct}, 
               the recapture probabilities for each capture occasion.} 
}
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\seealso{ \code{\link{closedp}}, \code{\link{closedpCI.t}}}
\examples{
data(hare)
closedp.Mtb(hare)

## Example producing an unstable estimate
data(mvole)
period4<-mvole[,16:20]
closedp.Mtb(period4)
}
\keyword{models}
