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
  \item{X}{ The table of the observed capture histories in one of the two accepted formats. In the default format, it has one row per unit captured in the experiment. In this case, the number of columns in the table represents the number of capture occasions in the experiment (noted \eqn{t}). In the alternative format, it contains one row per capture history followed by its frequency. In that case, \code{X} has \eqn{t}+1 columns. The first \eqn{t} columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. }
  \item{dfreq}{ This argument specifies the format of the data matrix \code{X}. By default, it is set to FALSE, which means that \code{X} has one row per unit. If it is set to TRUE, then the matrix \code{X} contains frequencies in its last column.}
  \item{x}{ An object, produced by the \code{closedp.Mtb} function, to print.}
  \item{\dots}{ Further arguments passed to or from other methods.}
  }
\details{ 
  The Mtb model is non-linear. It is fitted with the \code{optim} function instead of the \code{glm} fonction. Therefore, the abundance estimate can be unstable. 

  For the model to be identifiable, the parameters are constrained in the following way: \eqn{logit(c_i)=logit(p_i)+b}{logit(ci)=logit(pi)+b} for i in \eqn{2,\ldots,l}. 
}
\value{
\item{n}{ The number of captured units}
\item{results}{ A table containing the estimated population size, the standard error of estimation, the deviance, the number of degrees of freedom and the Akaike criteria.}
\item{parMtb}{ Capture-recapture parameters estimates for model Mtb : the abundance N, \eqn{p_1}{p1} to \eqn{p_t}{pt}, the probabilities of first capture for each capture occasion, and \eqn{c_2}{c2} to \eqn{c_t}{ct}, the recapture probabilities for each capture occasion.} 
}
\references{
Baillargeon, S. and Rivest, L.P. (2007). Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), http://www.jstatsoft.org/
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\seealso{ \code{\link{closedp}}, \code{\link{closedp.mX}}, \code{\link{closedp.h}} }
\examples{
data(hare)
closedp.Mtb(hare)

## Example producing an unstable estimate
data(mvole)
period4<-mvole[,16:20]
closedp.Mtb(period4)
}
\keyword{models}
