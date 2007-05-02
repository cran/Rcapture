\name{robustd}
\alias{robustd.t}
\alias{robustd.0}
\alias{print.robustd}
\title{ Loglinear Models for Capture-Recapture Experiments Carried Out According to a Robust Design}
\description{
These functions compute various demographic parameters and capture probabilities per period using loglinear robust design models in capture-recapture experiments.
\code{robustd.t} and \code{robustd.0} fit the model on different response variable. 
\code{robustd.t} uses the frequencies of the observable capture histories in terms of capture success or failure for each capture occasions of each primary period (generated by \code{histfreq.t}).
\code{robustd.0} uses the frequencies of the observable capture histories in terms of number of captures per primary period (generated by \code{histfreq.0}).
}
\usage{
robustd.t(X, dfreq = FALSE, vt, vm = "M0", vh = list("Chao"),
          va = 2, neg = TRUE)

robustd.0(X, dfreq = FALSE, vt, vm = "M0", vh = list("Chao"), 
          va = 2, neg = TRUE)

\method{print}{robustd}(x, \dots)
}
\arguments{
  \item{X}{ The table of the observed capture histories in one of the two accepted formats. In the default format, it has one row per unit captured in the experiment and sum(vt) columns. In the alternative format, it contains one row per capture history followed by its frequency. In that case, \code{X} has sum(vt)+1 columns. The first sum(vt) columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. }
  \item{dfreq}{ This argument specifies the format of the data matrix \code{X}. By default, it is set to FALSE, which means that \code{X} has one row per unit. If it is set to TRUE, then the matrix \code{X} contains frequencies in its last column.}
  \item{vt}{ A vector containing the numbers of capture occasions for each primary sampling period. The length of this vector equals the number of primary sampling periods (noted \eqn{I}). }
  \item{vm}{ A vector indicating the closed population model for each primary period. The elememts of \code{vm} can be "none"=no model, "M0"=M0 model, "Mt"=Mt model, "Mh"=Mh model or "Mth"=Mth model for \code{robustd.t}. For \code{robustd.0}, models with time effect are not allowed. So \code{vm} can than be "none"=no model, "M0"=M0 model, or "Mh"=Mh model. 
 The 'no model' cannot be selected for the first or the last period. If a single character string is given for \code{vm}, the corresponding model is used for all periods. The default value for \code{vm} is "M0" for all period.}
  \item{vh}{ A list indicating, for each primary period with a heterogeneity model, the form of the columns for heterogeneity in the design matrix. The elements of \code{vht} can be "Chao", "Poisson", "Darroch" or any numerical \code{R} function created beforehand by the user. "Chao" represents Chao's model, "Poisson" represents the function \eqn{f(k)=a^k-1}, where \eqn{k} is the number of captures, and "Darroch" represents the function \eqn{f(k)=k^2/2}. If an \code{R} function is given, it is the implemantation of any mathematical function \eqn{f(k)}. It has only one argument. For the Poisson model, the exponent's base \code{a} is specified in the argument \code{va}. }
  \item{va}{ A vector indicating, for each primary period with a Poisson heterogeneity model, the value of the exponent's base in \eqn{f(k)=a^k-1}. If  \code{va} consists of a single number, this number is used for all the periods with a Poisson model. The default value for \code{va} is 2 for all the periods with a Poisson model. }
  \item{neg}{ If this option is set to TRUE, negative gamma parameters and negative eta parameters in Chao's models are set to zero. This insures that the estimated survival probabilities belong to [0, 1] and that the births are positive.}
  \item{x}{ An object, produced by the \code{robustd.t} or the \code{robustd.0} function, to print.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{
These functions also generate statistics to test the presence of temporary emigration.

The Poisson regression used to fit a robust design model has one entry for each possible caputre history, including those that are unobserved.  The size of the dependent vector is therefore \eqn{2^{\mbox{sum(vn)}}-1}{2^sum(vn)-1} for \code{robustd.t}. Models with a large sum(vt) are hard to fit with \code{robustd.t}. \code{robustd.0} uses a more parsimonious coding for the capture histories and can fit larger models.

Standard errors are calculated by linearization.
}
\value{
  \item{n}{ The number of captured units}
  \item{models }{ A vector of length \eqn{I} identifying the closed population models chosen for each period.}
  \item{model.fit }{ A table containing the deviance, degrees of freedom and AIC of the fitted model.}
  \item{emig.fit }{ A table containing, for the model with an added temporary emigration effect, the deviance, degrees of freedom and AIC.}
  \item{emig.param }{ The estimated temporary emigration parameters and their standards errors. The \eqn{I-1} first rows are estimations of the differences \eqn{logit(p^b)-logit(p^w)} for periods 2 to \eqn{I-1} (\eqn{p^b} represents a between primary period and \eqn{p^w} a within primary period estimate of the capture probability). The last row gives a pooled estimate of these differences calculated under the assumption that they are homogenous. Negative estimates are associated with a temporary emigration.}
  \item{capture.prob }{ The estimated capture probabilities per period and their standard errors.}
  \item{survivals }{ The estimated survival probabilities between periods and their standard errors.}
  \item{N }{ The estimated population sizes per period and their standard errors. }
  \item{birth}{ The estimated number of new arrivals in the population between periods and their standard errors.}
  \item{Ntot}{ The estimated total number of units who ever inhabited the survey area and its standard error.} 
  \item{loglin.param }{The loglinear model parameters estimations and their standard errors, calculated by the \code{glm} function. }
  \item{u.vector }{The Ui statistics, useful for the survival probabilities calculation, and their standard errors}
  \item{v.vector }{The Vi statistics, useful for the population sizes estimation, and their standard errors}
  \item{cov}{ The covariance matrix of all the demographic parameters estimates. }
  \item{neg}{ The position of the gamma and eta parameters set to zero in the loglinear parameter vector.}
}
\references{ 
Baillargeon, S. and Rivest, L.P. (2007). The Rcapture package: Loglinear models for capture-recapture in R. \emph{Journal of Statistical Software}, to appear (available online at \url{http://www.jstatsoft.org/}).

Rivest, L.P. and Daigle, G. (2004) Loglinear models for the robust design in mark-recapture experiments. \emph{Biometrics}, \bold{60}, 100--107.
}
\author{Sophie Baillargeon \email{sbaillar@mat.ulaval.ca} and Louis-Paul Rivest \email{lpr@mat.ulaval.ca} }
\note{ 
This function uses the \code{glm} function of the \code{stats} package. 
}
\seealso{\code{\link{closedp}}, \code{\link{openp}} }
\examples{
data(mvole)

  # The mvole data set contains a total of 30 capture occasions (the 
  # tenth capture occasion doesn't have any new capture and is taken
  # out of the analysis). This number being large, we can only use 
  # the robustd.0 function to fit a robust design model.
robustd.0(mvole[,-10],vt=c(5,4,rep(5,4)),vm="Mh",vh="Poisson",va=1.5)
# Not run: 
# robustd.t(mvole[,-10],vt=c(5,4,rep(5,4)),vm="Mh",vh="Poisson",va=1.5)
# should fail

  # Considering only the first 3 periods of the data set, we can use the 
  # robustd.t function to fit a model with a temporal effect.
robustd.t(mvole[,c(1:9,11:15)],vt=c(5,4,5),vm="Mth",vh="Poisson",va=1.5)
}
\keyword{ models }