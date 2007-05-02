\name{openp}
\alias{openp}
\alias{print.openp}
\alias{plot.openp}
\title{ Loglinear Models for Open Population Capture-Recapture Experiments  }
\description{
  This function computes various demographic parameters using a loglinear model for open populations in capture-recapture experiments.
}
\usage{
openp(X, dfreq=FALSE, m="up" , neg=TRUE, keep=rep(TRUE,2^I-1))

\method{print}{openp}(x, \dots)

\method{plot}{openp}(x, \dots)
}
\arguments{
  \item{X}{ The table of the observed capture histories in one of the two accepted formats. In the default format, it has one row per unit captured in the experiment. In this case, the number of columns in the table represents the number of capture occasions in the experiment (noted \eqn{I}). In the alternative format, it contains one row per capture history followed by its frequency. In that case, \code{X} has \eqn{I}+1 columns. The first \eqn{I} columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. }
  \item{dfreq}{ This argument specifies the format of the data matrix \code{X}. By default, it is set to FALSE, which means that \code{X} has one row per unit. If it is set to TRUE, then the matrix \code{X} contains frequencies in its last column.}
  \item{m}{ This argument is a character string taking the value "up" (up = unconstrained probabilities) or "ep" (ep = equal probabilities). If \code{m} is set to "up" (the default), no constraint is fixed on the loglinear parameters. Therefore some loglinear and demographic parameters are not estimable. On the other hand, when \code{m} takes the value "ep", the capture probabilities are set to a common value; this enables the estimation of all the parameters.}
  \item{keep}{ This option is useful to fit the model on a subset of the possible capture histories. \code{keep} is a logical vector of length \eqn{2^I-1} taking the value TRUE for a history kept and FALSE for a history put aside. In this vector, the order of the capture histories is as defined in the \code{histpos.t} function. By default, every capture history is kept. }
  \item{neg}{ If this option is set to TRUE, relevant negative gamma parameters are set to zero. This insures that the estimated survival probabilities belong to [0, 1] and that the births are positive.}
  \item{x}{ An object, produced by the \code{openp} function, to print or to plot.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
The function \code{openp} generates statistics to test the presence of a trap effect.

The \code{plot.openp} function produces a scatterplot of the Pearson residuals of the model versus the frequencies of capture.

If the data matrix \code{X} was obtained through the \code{periodhist} function, the \code{dfreq} argument must be set to TRUE.

Standard errors are calculated by linearization.
}
\value{
  \item{n}{ The number of captured units}
  \item{model.fit }{ A table containing the deviance, degrees of freedom and AIC of the fitted model.}
  \item{trap.fit }{ A table containing, for the models with an added trap effect, the deviance, degrees of freedom and AIC.}
  \item{trap.param }{ The estimated trap effect parameters and their standard errors. For m="up", the \eqn{I-3} first rows of \code{trap.param} are estimations of the differences \eqn{logit}(capture probability after a capture)-\eqn{logit}(capture probability after a miss) for periods 3 to \eqn{I-1}. The last row gives a pooled estimate of these differences calculated under the assumption that they are homogenous.}
  \item{capture.prob }{ The estimated capture probabilities per period and their standard errors.}
  \item{survivals }{ The estimated survival probabilities between periods and their standard errors.}
  \item{N }{ The estimated population sizes per period and their standard errors. }
  \item{birth}{ The estimated number of new arrivals in the population between periods and their standard errors.}
  \item{Ntot}{ The estimated total number of units who ever inhabited the survey area and its standard error.} 
  \item{glm}{ The 'glm' object obtained from fitting the loglinear model}
  \item{loglin.param }{The loglinear model parameters estimations and their standard errors, calculated by the \code{glm} function. }
  \item{u.vector }{The Ui statistics, useful for the survival probabilities calculation, and their standard errors}
  \item{v.vector }{The Vi statistics, useful for the population sizes estimation, and their standard errors}
  \item{cov}{ The covariance matrix of all the demographic parameters estimates. }
  \item{neg}{ The position of the gamma parameters set to zero in the loglinear parameter vector.}
}
\references{ 
Baillargeon, S. and Rivest, L.P. (2007). The Rcapture package: Loglinear models for capture-recapture in R. \emph{Journal of Statistical Software}, to appear (available online at \url{http://www.jstatsoft.org/}).

Rivest, L.P. and Daigle, G. (2004) Loglinear models for the robust design in mark-recapture experiments. \emph{Biometrics}, \bold{60}, 100--107.
}
\author{Sophie Baillargeon \email{sbaillar@mat.ulaval.ca} and Louis-Paul Rivest \email{lpr@mat.ulaval.ca} }
\note{ 
If your data contains more than one capture occasion within primary periods, use the \code{\link{periodhist}} function to create the input data matrix \code{X} needed by the \code{openp} function.

This function uses the \code{glm} function of the \code{stats} package. 
}
\seealso{\code{\link{closedp}}, \code{\link{periodhist}}, \code{\link{robustd.t}}, \code{\link{robustd.0}}}
\examples{
data(duck)
op.m1 <- openp(duck, dfreq=TRUE)
plot(op.m1)

# To remove the capture history 111111.
keep2 <- apply(histpos.t(6),1,sum)!=6
op.m2 <- openp(duck, dfreq=TRUE, keep=keep2)
op.m2

# To remove the capture histories with 5 captures or more
keep3 <- apply(histpos.t(6),1,sum)<5
op.m3 <- openp(duck, dfreq=TRUE, keep=keep3)
op.m3


data(mvole)
mvole.op<-periodhist(mvole,vt=rep(5,6))
openp(mvole.op, dfreq=TRUE)
}
\keyword{ models }
