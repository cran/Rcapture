\name{closedp.custom}
\alias{closedp.custom}
\alias{closedp.mX}
\alias{closedp.h}
\alias{print.closedp.custom}
\alias{boxplot.closedp.custom}
\title{ Customization of a Loglinear Model for Closed Populations Capture-Recapture Experiments }
\description{
This function is deprecated, please use \code{closedpCI.t} instead. 

The \code{closedp.mX} function fits a loglinear model given a design matrix \code{mX}. The \code{closedp.h} function fits Mh or Mth models for which the form of the column for heterogeneity in the design matrix is determined by the user. 
}
\usage{
closedp.mX(X, dfreq=FALSE, mX, mname="Customized model")

closedp.h(X, dfreq=FALSE, m="Mh", h="Poisson", a=2)

\method{print}{closedp.custom}(x, \dots)

\method{boxplot}{closedp.custom}(x, \dots)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} for a description of the accepted formats).  }
  \item{dfreq}{ A logical. By default FALSE, which means that \code{X} has one row per unit. If TRUE, it indicates that the matrix \code{X} contains frequencies in its last column.}
  \item{mX}{ The design matrix of the loglinear model. In this matrix, the order of the capture histories is as defined in the \code{histpos.t} function.}
  \item{mname}{ A character string specifying the name of the customized model.}
  \item{m}{ A character string indicating the model to fit, either "Mh"=Mh model or "Mth"=Mth model}
  \item{h}{ The character string "Poisson" (\eqn{2^k-1}) or a numerical function taking only one argument specifying the form of the column for heterogeneity in the design matrix.}
  \item{a}{ The value of the exponent's base for a Poisson model. }
  \item{x}{ An object, produced by the \code{closedp.mX} or the \code{closedp.h} function, to print or to plot.}
  \item{\dots}{ Further arguments passed to or from other methods.}
 }
\details{ 
  An intercept is added to the model. Therefore, the mX matrix must not contain a column of ones.
  
  The abundance estimation is calculated as the number of captured units plus the exponential of the intercept. Therefore, these functions are not suited for models with a behavioral effect.

  In \code{closedp.h}, the argument \code{h} cannot take the value "Chao" or "Darroch" because these models are already fitted by the \code{closedp} function.

  The \code{boxplot.closedp.custom} function produces a boxplot of the pearson residuals of the customized model.
}
\value{
\item{n}{ The number of captured units}
\item{results}{ A table containing the estimated population size, the standard error of estimation, the deviance, the number of degrees of freedom and the Akaike criteria.}
\item{glm}{ The 'glm' object obtained from fitting the model.} 
}
\references{ 
Rivest, L.P. and Baillargeon, S. (2007) Applications and extensions of Chao's moment estimator for the size of a closed population. \emph{Biometrics}, \bold{63}(4), 999--1006.
}
\author{Sophie Baillargeon \email{sbaillar@mat.ulaval.ca} and Louis-Paul Rivest \email{lpr@mat.ulaval.ca} }
\note{ These functions use the \code{glm} function of the \code{stats} package. }
\seealso{\code{\link{closedp}}, \code{\link{closedp.Mtb}} }
\examples{
data(HIV)
mat<-histpos.t(4)
mX2<-cbind(mat,mat[,1]*mat[,2])
closedp.mX(HIV,dfreq=TRUE,mX=mX2)

data(mvole)
period3<-mvole[,11:15]
psi <- function(x) { -log(3.5 + x) + log(3.5) }
closedp.h(period3, h = psi)
}
\keyword{ models }
