\name{closedp.Mtb}
\alias{closedp.Mtb}
\alias{print.closedp.Mtb}
\title{ Closed Population Capture-Recapture Model Mtb }
\description{
  This function fits model Mtb for closed populations in capture-recapture experiments.
}
\usage{
closedp.Mtb(X, dfreq=FALSE, method = "BFGS", \dots)

\method{print}{closedp.Mtb}(x, \dots)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} for a description of the accepted formats).  }
  \item{dfreq}{ A logical. By default \code{FALSE}, which means that \code{X} has one row per unit. If TRUE, it 
                indicates that the matrix \code{X} contains frequencies in its last column.}
  \item{method}{ The method to be used by \code{\link{optim}}. The default is \code{"BFGS"}.}
  \item{\dots}{ Further arguments to be passed to \code{\link{optim}} or \code{\link{print.default}}.}
  \item{x}{ An object, produced by the \code{closedp.Mtb} function, to print.}
  }
\details{ 
  The Mtb model is non-linear. It is fitted with the \code{\link{optim}} function instead of the \code{glm} function. Therefore, the abundance estimate can be unstable. 

  For the model to be identifiable, the parameters are constrained in the following way: 
  \eqn{logit(c_i)=logit(p_i)+b}{logit(ci)=logit(pi)+b} for i in \eqn{2,\ldots,l}. 
}
\value{
\item{n}{ The number of captured units}
\item{t}{ The total number of capture occasions in the data matrix \code{X}.}
\item{results}{ A table containing, for the fitted model:
    \describe{
    \item{\code{abundance}:}{  the estimated population size,}
    \item{\code{stderr}:}{  the standard error of the estimated population size,}
    \item{\code{deviance}:}{  the model's deviance,}
    \item{\code{df}:}{  the number of degrees of freedom,}
    \item{\code{AIC}:}{  the Akaike's information criterion,}
    \item{\code{BIC}:}{  the bayesian information criterion,}
    \item{\code{infoFit}:}{  a numerical code giving information about error or warnings 
          encountered when fitting the model (see \code{\link{Rcapture-package}} for details).}
    }}
\item{optim}{The output produced by \code{\link{optim}} from fitting the model.}
\item{optim.warn}{ A vector of character strings. If the \code{\link{optim}} function generates
                   one or more warnings when fitting the model, a copy of these warnings are
                   stored in \code{optim.warn}. \code{NULL} if \code{\link{optim}} did not produce
                   any warnings. }
\item{parMtb}{ Capture-recapture parameters estimates for model Mtb : the abundance N, \eqn{p_1}{p1} to \eqn{p_t}{pt},
               the probabilities of first capture for each capture occasion, and \eqn{c_2}{c2} to \eqn{c_t}{ct}, 
               the recapture probabilities for each capture occasion.} 
}
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \doi{10.18637/jss.v019.i05}.
}
\author{ Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and Sophie Baillargeon }

\seealso{ \code{\link{closedp}}, \code{\link{closedpCI.t}}}
\examples{
# hare data set
closedp.Mtb(hare)

## Example producing an unstable estimate
# Fourth primary period of mvole data set
period4 <- mvole[, 16:20]
closedp.Mtb(period4)
}
\keyword{models}
