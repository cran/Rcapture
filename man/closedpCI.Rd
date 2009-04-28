\name{closedpCI}
\alias{closedpCI}
\alias{closedpCI.t}
\alias{closedpCI.0}
\alias{print.closedpCI}
\alias{plotCI}
\alias{plotCI.closedpCI}
\alias{boxplot.closedpCI}
\alias{plot.closedpCI}
\title{ Customization of a Loglinear Model and profile Likelihood Confidence Interval for Abundance Estimation in Closed Population Capture-Recapture Experiments }
\description{
   The \code{closedpCI.t} and \code{closedpCI.0} functions fit a loglinear model specified by the user and computes the multinomial profile likelihood confidence interval for the adundance estimation. The model can be given as a design matrix \code{mX} or identified trougth arguments \code{m}, \code{h} and \code{theta}. These functions extand \code{closedp.t} and \code{closedp.0} as they broaden the range of model one can fit and they computes confidence interval. Unlike the \code{closedp} functions, it fits only one model at a time. 
}
\usage{
closedpCI.t(X, dfreq=FALSE, m=c("M0","Mt","Mh","Mth"), 
            h=c("Chao","Poisson","Darroch","Gamma"), theta=2, 
            mX=NULL, mname, neg=TRUE, alpha=0.05)

closedpCI.0(X, dfreq=FALSE, dtype=c("hist","nbcap"), t, t0=t, 
            m=c("M0","Mh"), h=c("Chao","Poisson","Darroch","Gamma"), 
            theta=2, mX=NULL, mname, neg=TRUE, alpha=0.05)
			
\method{print}{closedpCI}(x, \dots)

plotCI(x, \dots)
\method{plotCI}{closedpCI}(x, main="Profile Likelihood Confidence Interval", \dots)

\method{boxplot}{closedpCI}(x, main="Boxplots of Pearson Residuals", \dots)

\method{plot}{closedpCI}(x, main="Scatterplot of Pearson Residuals", \dots)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} for a description of the accepted formats).  }
  \item{dfreq}{ A logical. By default FALSE, which means that \code{X} has one row per unit. If TRUE, it indicates that the matrix \code{X} contains frequencies in its last column.}
  \item{dtype}{ A characters string, either "hist" or "nbcap", to specify the type of data. "hist", the default, means that \code{X} contains complete observed capture histories. "nbcap" means that \code{X} contains numbers of captures (see \code{\link{Rcapture-package}} for details on data formats). }
  \item{t}{ Requested only if \code{dtype="nbcap"}. A numeric specifying the total number of capture occasions in the experiment.}
  \item{t0}{ A numeric. Models are fitted considering only the frequencies of units captured 1 to \code{t0} times. By default \code{t0=t}. }
  \item{m}{ A character string indicating the model to fit. For \code{closedpCI.0} it can be either "M0"=M0 model or "Mh"=Mh model. For \code{closedpCI.t} it can also be "Mt"=Mt model or "Mth"=Mth model.}
  \item{h}{ A character string ("Chao", "Poisson", "Darroch" or "Gamma") or a numerical \code{R} function specifying the form of the column for heterogeneity in the design matrix. "Chao" represents Chao's model, "Poisson" represents the function \eqn{f(k)=theta^k-1}, where \eqn{k} is the number of captures, "Darroch" represents the function \eqn{f(k)=k^2/2}, and "Gamma" represents the function \eqn{f(k)=-log(theta + k) + log(theta)}. If an \code{R} function is given, it is the implemantation of any convex mathematical function \eqn{f(k)}. It has only one argument.}
  \item{theta}{ The value of the parameter for a Poisson or Gamma model.}
  \item{mX}{  The design matrix of the loglinear model. In this matrix, the order of the capture histories is as defined in the \code{histpos.t} or \code{histpos.0} function.}
  \item{mname}{ A character string specifying the name of the customized model.}
  \item{neg}{ If this option is set to TRUE, negative eta parameters in Chao's models are set to zero. }
  \item{alpha}{ A confidence interval with confidence level 1-alpha is constructed. The value of alpha must be between 0 and 1; the default is 0.05.  }
  \item{x}{ An object, produced by the \code{closedpCI.t} function, to print.}
  \item{main}{ A main title for the plot}
  \item{\dots}{ Further arguments to be passed to methods (see \code{print.default}, \code{plot.default} or \code{boxplot.default}).}
  }
\details{ 
The \code{closedpCI.t} function fits models using the frequencies of the observable capture histories (vector of size \eqn{2^t-1}), whereas \code{closedp.0} uses the number of units capture i times, for \eqn{i=1,\ldots,t}{i=1,\ldots,t} (vector of size \eqn{t}). Thus, \code{closedpCI.0} can be used with data sets larger than those for \code{closedpCI.t}. 

   This function does not work for closed population models featuring a behavioral effect, such as Mb and Mbh. The abundance estimation is calculated as the number of captured units plus the exponential of the Poisson regression intercept. However, models with a behavioral effect can by fitted with \code{\link{closedp.t}} (Mb and Mbh), \code{\link{closedp.Mtb}} and \code{\link{closedp.bc}}.
  
  An intercept is added to the model. Therefore, the \code{mX} matrix must not contain a column of ones.

  The \code{plotCI.closedpCI} function produces a plot of the multinomial profile likelihood for N. The value of N maximizing the profile likelihood and the bounds of the confidence interval are identified.
  
  The \code{boxplot.closedpCI} function produces a boxplot of the Pearson residuals of the customized model.

  The \code{plot.closedpCI} function traces the scatterplot of the Pearson residuals in terms of \eqn{f_i}{fi} (number of units captured i times) for the customized model.  
  }
\value{
\item{n}{ The number of captured units}
\item{t}{ The number of capture occasions in the data matrix \code{X}}
\item{results}{ A table containing the estimated population size, the standard error of estimation, the deviance, the number of degrees of freedom and the Akaike criteria.}
\item{glm}{ The 'glm' object obtained from fitting the model.} 
\item{CI}{ A table containing the abundance estimation and its confidence interval.}
\item{alpha}{ 1-the confidence level of the interval.}
\item{NCI}{ The x-coordinates for \code{plot.closedpCI.t}}
\item{loglikCI}{ The y-coordinates for \code{plot.closedpCI.t}} 
\item{t0}{ A copy of the \code{t0} argument given in the function call.}
}
\references{ 
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.

Rivest, L.P. and Baillargeon, S. (2007) Applications and extensions of Chao's moment estimator for the size of a closed population. \emph{Biometrics}, \bold{63}(4), 999--1006.

Cormack, R. M. (1992) Interval estimation for mark-recapture studies of closed populations. \emph{Biometrics}, \bold{48}, 567--576.
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\note{ This function uses the \code{optimize} and the \code{uniroot} functions of the \code{stats} package. }
\seealso{\code{\link{closedp}},  \code{\link{closedp.Mtb}}}
\examples{
data(hare)
CI<-closedpCI.t(hare, m = "Mth", h = "Poisson", theta = 2)
CI
plotCI(CI)

data(HIV)
mat<-histpos.t(4)
mX2<-cbind(mat,mat[,1]*mat[,2])
closedpCI.t(HIV,dfreq=TRUE,mX=mX2,mname="Mt interaction 1,2")

data(BBS2001)
CI0<-closedpCI.0(BBS2001,dfreq=TRUE,dtype="nbcap",t=50,t0=20,
                 m="Mh",h="Gamma",theta=3.5)
CI0
plot(CI0)
plotCI(CI0)

### As an alternative to a gamma model, one can fit a negative Poisson model.
### It is appropriate in experiments where very small capture probabilities
### are likely. It can lead to very large estimators of abundance. 
data(mvole)
period3<-mvole[,11:15]
psi <- function(x) { 0.5^x - 1 }
closedpCI.t(period3, m = "Mh", h = psi)
}
\keyword{ models }
