\name{profileCI}
\alias{profileCI}
\alias{print.profileCI}
\title{ Profile Likelihood Confidence Interval for Abundance Estimation in Closed Populations Capture-Recapture Experiments }
\description{
This function is deprecated, please use \code{closedpCI.t} instead.

The \code{profileCI} function computes the multinomial profile likelihood for the abundance of some closed population capture-recapture models. 
}
\usage{
profileCI(X, dfreq=FALSE, m="M0", h="Chao", a=2, mX=NULL, 
          mname="Customized model", neg=TRUE, alpha=0.05)

\method{print}{profileCI}(x, \dots)
}
\arguments{
  \item{X}{ The table of the observed capture histories in one of the two accepted formats. In the default format, it has one row per unit captured in the experiment. In this case, the number of columns in the table represents the number of capture occasions in the experiment (noted \eqn{t}). In the alternative format, it contains one row per capture history followed by its frequency. In that case, \code{X} has \eqn{t}+1 columns. The first \eqn{t} columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. }
  \item{dfreq}{ This argument specifies the format of the data matrix \code{X}. By default, it is set to "FALSE", which means that \code{X} has one row per unit. If it is set to TRUE, then the matrix \code{X} contains frequencies in its last column.}
  \item{m}{ A character string identifying the  model, either "M0"=M0 model, "Mt"=Mt model, "Mh"=Mh model or "Mth"=Mth model.}
  \item{h}{ A character string ("Chao" = Chao model, "Poisson" = \eqn{2^k-1}), "Darroch" = \eqn{(k^2)/2} or a numerical function taking only one argument specifying the form of the column for heterogeneity in the design matrix.}
  \item{a}{ The value of the exponent's base for a Poisson model. }
  \item{mX}{  The design matrix of the loglinear model. In this matrix, the order of the capture histories is as defined in the \code{histpos.t} function.}
  \item{mname}{ A character string specifying the name of the customized model.}
  \item{neg}{ If this option is set to TRUE, negative eta parameters in Chao models are set to zero. }
  \item{alpha}{ A confidence interval with confidence level 1-alpha is constructed. The value of alpha must be between 0 and 1; the default is 0.05.  }
  \item{x}{ An object, produced by the \code{profileCI} function, to print.}
  \item{\dots}{ Further arguments passed to or from other methods.}
  }
\details{ 
  This function does not work for closed population models featuring a behavioral effect, such as Mb and Mbh.
}
\value{
This function produces a plot of the multinomial profile likelihood for N. The value of N maximizing the profile likelihood and the bounds of the confidence interval are identified.
It also produces the following objects :
\item{n}{ The number of captured units}
\item{results}{ A table containing the abundance estimation and its confidence interval.}
\item{alpha}{ 1-the confidence level of the interval.}
}
\references{ 
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.

Cormack, R. M. (1992) Interval estimation for mark-recapture studies of closed populations. \emph{Biometrics}, \bold{48}, 567--576.
 }
\author{Sophie Baillargeon \email{sbaillar@mat.ulaval.ca} and Louis-Paul Rivest \email{lpr@mat.ulaval.ca} }
\seealso{\code{\link{closedp}}, \code{\link{closedp.mX}}, \code{\link{closedp.h}} }
\examples{
data(hare)
profileCI(hare, m = "Mth", h = "Poisson", a = 2)

data(HIV)
mat<-histpos.t(4)
mX2<-cbind(mat,mat[,1]*mat[,2])

profileCI(hare,m="Mh",h="Chao")
}
\keyword{ models }
