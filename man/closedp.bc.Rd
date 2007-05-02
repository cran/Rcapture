\name{closedp.bc}
\alias{closedp.bc}
\alias{print.closedp.bc}
\title{ Bias Correction for Closed Population Loglinear Models }
\description{
  This function applies a bias correction to the abundance estimations obtained by the \code{closedp} function.
}
\usage{
closedp.bc(X, dfreq=FALSE)

\method{print}{closedp.bc}(x, \dots)
}
\arguments{
  \item{X}{ The table of the observed capture histories in one of the two accepted formats. In the default format, it has one row per unit captured in the experiment. In this case, the number of columns in the table represents the number of capture occasions in the experiment (noted \eqn{t}). In the alternative format, it contains one row per capture history followed by its frequency. In that case, \code{X} has \eqn{t}+1 columns. The first \eqn{t} columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. }
  \item{dfreq}{ This argument specifies the format of the data matrix \code{X}. By default, it is set to FALSE, which means that \code{X} has one row per unit. If it is set to TRUE, then the matrix \code{X} contains frequencies in its last column.}
  \item{x}{ An object, produced by the \code{closedp.bc} function, to print.}
  \item{\dots}{ Further arguments passed to or from other methods.}
  }
\details{ 
The bias correction is done through frequency modifications as described in Rivest and Levesque (2001). The variances calculated with the modified frequencies are less biased than the standard ones, but they can overestimate the mean squared errors, especially when the data is sparse.
}
\value{
\item{n}{ The number of captured units}
\item{results}{ A table containing, for each model, the corrected population size estimation and its standard error.}
}
\references{ 
Baillargeon, S. and Rivest, L.P. (2007). The Rcapture package: Loglinear models for capture-recapture in R. \emph{Journal of Statistical Software}, to appear (available online at \url{http://www.jstatsoft.org/}).

Rivest, L.P. and Levesque, T. (2001). Improved log-linear model estimators of abundance in capture-recapture experiments. \emph{Canadian Journal of Statistics}, \bold{29}, 555--572.
}
\author{Sophie Baillargeon \email{sbaillar@mat.ulaval.ca} and Louis-Paul Rivest \email{lpr@mat.ulaval.ca} }
\note{ This function uses the \code{glm} function of the \code{stats} package. }
\seealso{\code{\link{closedp}} }
\examples{
data(mvole)
period3<-mvole[,11:15]
closedp.bc(period3)
}
\keyword{ models }
