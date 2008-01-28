\name{closedp}
\alias{closedp}
\alias{print.closedp}
\alias{boxplot.closedp}
\title{ Loglinear Models for Closed Population Capture-Recapture Experiments }
\description{
  This function fits various loglinear models for closed populations in capture-recapture experiments: M0, Mt, Mh Chao, Mh Poisson2, Mh Darroch, Mth Chao, Mth Poisson2, Mth Darroch, Mb and Mbh.
}
\usage{
closedp(X, dfreq=FALSE, neg=TRUE)

\method{print}{closedp}(x, \dots)

\method{boxplot}{closedp}(x, \dots)
}
\arguments{
  \item{X}{ The table of the observed capture histories in one of the two accepted formats. In the default format, it has one row per unit captured in the experiment. In this case, the number of columns in the table represents the number of capture occasions in the experiment (noted \eqn{t}). In the alternative format, it contains one row per capture history followed by its frequency. In that case, \code{X} has \eqn{t}+1 columns. The first \eqn{t} columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. }
  \item{dfreq}{ This argument specifies the format of the data matrix \code{X}. By default, it is set to FALSE, which means that \code{X} has one row per unit. If it is set to TRUE, then the matrix \code{X} contains frequencies in its last column.}
  \item{neg}{ If this option is set to TRUE, negative eta parameters in Chao's models are set to zero. }
  \item{x}{ An object, produced by the \code{closedp} function, to print or to plot.}
  \item{\dots}{ Further arguments passed to or from other methods.}
  }
\details{ 
  Chao's models estimate a lower bound for the abundance, both with (Mth Chao) and without a time effect (Mh Chao). The estimate obtained under Mh Chao is Chao's (1987) moment estimator. Rivest and Baillargeon (2007) exhibit a loglinear model underlying this estimator and provide a generalization to Mth. For these two models, a small deviance means that there is an heterogeneity in capture probabilities; it does not mean that the lower bound estimates are unbiased.
For Darroch's models, the column for heterogeneity in the design matrix is defined as \eqn{k^2/2} where \eqn{k} is the number of captures;  these models for Mh and Mth are considered by Darroch et al. (1993) and Agresti (1994). For Poisson2 models, the column for heterogeneity in the design matrix is \eqn{2^k-1}. Poisson models with an exponent's base different than 2 can be fitted with the \code{closedp.h} function. These models are discussed in Rivest and Baillargeon (2007);  they typically yield smaller corrections for heterogeneity than Darroch's model since the capture probabilities are bounded from below under these models.
  
  The \code{boxplot.closedp} function produces boxplots of the Pearson residuals of the ten loglinear models.

  To perform a bias correction on the abundance estimates, use the \code{closedp.bc} function.
}
\value{
\item{n}{ The number of captured units}
\item{t}{ The number of capture occasions in the data matrix \code{X}}
\item{results}{ A table containing, for each model, the estimated population size, the standard error of estimation, the deviance, the number of degrees of freedom and the Akaike criteria.}
\item{glmM0}{ The 'glm' object obtained from fitting the M0 model.} 
\item{glmMt}{ The 'glm' object obtained from fitting the Mt model.} 
\item{glmMhC}{ The 'glm' object obtained from fitting the Mh Chao model.} 
\item{glmMhP}{ The 'glm' object obtained from fitting the Mh Poisson model.} 
\item{glmMhD}{ The 'glm' object obtained from fitting the Mh Darroch model.} 
\item{glmMthC}{ The 'glm' object obtained from fitting the Mth Chao model.} 
\item{glmMthP}{ The 'glm' object obtained from fitting the Mth Poisson model.} 
\item{glmMthD}{ The 'glm' object obtained from fitting the Mth Darroch model.} 
\item{glmMb}{ The 'glm' object obtained from fitting the Mb model.} 
\item{glmMbh}{ The 'glm' object obtained from fitting the Mbh model.}
\item{parM0}{ Capture-recapture parameters estimates for model M0 : the abundance N and p, the capture probability at any capture occasion.} 
\item{parMt}{ Capture-recapture parameters estimates for model Mt : N and \eqn{p_1}{p1} to \eqn{p_t}{pt}, the capture probabilities for each capture occasion.} 
\item{parMhC}{ Capture-recapture parameters estimates for model MhC : N and p, the average probability of capture.} 
\item{parMhP}{ Capture-recapture parameters estimates for model MhP : N and p, the average probability of capture.} 
\item{parMhD}{ Capture-recapture parameters estimates for model MhD : N and p, the average probability of capture.} 
\item{parMthC}{ Capture-recapture parameters estimates for model MthC : N and \eqn{p_1}{p1} to \eqn{p_t}{pt}, the average probabilities of capture for each occasion.} 
\item{parMthP}{ Capture-recapture parameters estimates for model MthP : N and \eqn{p_1}{p1} to \eqn{p_t}{pt}, the average probabilities of capture for each occasion.} 
\item{parMthD}{ Capture-recapture parameters estimates for model MthD : N and \eqn{p_1}{p1} to \eqn{p_t}{pt}, the average probabilities of capture for each occasion.} 
\item{parMb}{ Capture-recapture parameters estimates for model Mb : N, p, the probability of first capture at any capture occasion, and c, the recapture probability at any capture occation.} 
\item{parMbh}{ Capture-recapture parameters estimates for model Mbh : N, p and c.}
\item{negMhC}{ The position of the eta parameters set to zero in the loglinear parameter vector of model MhC.}
\item{negMthC}{ The position of the eta parameters set to zero in the loglinear parameter vector of model MthC.}
}
\references{ 
Agresti, A. (1994). Simple capture-recapture models permitting unequal catchability and variable sampling effort. \emph{Biometrics}, 50, 494--500.

Baillargeon, S. and Rivest, L.P. (2007). Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), http://www.jstatsoft.org/

Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchabililty.  \emph{ Biometrics}, 45, 427--438

Darroch, S.E.,  Fienberg, G.,  Glonek, B. and Junker, B. (1993). A three sample multiple capture-recapture approach to the census population estimation with heterogeneous catchability. \emph{Journal of the American Statistical Association}, 88, 1137--1148

Rivest, L.P. and Levesque, T. (2001) Improved log-linear model estimators of abundance in capture-recapture experiments. \emph{Canadian Journal of Statistics}, 29, 555--572.

Rivest, L.P. and Baillargeon, S. (2007). Applications and extensions of Chao's moment estimator for the size of a closed population. \emph{Biometrics}, \bold{63}(4), 999--1006.
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\note{ This function uses the \code{glm} function of the \code{stats} package. }
\seealso{\code{\link{uifit}}, \code{\link{closedp.bc}}, \code{\link{closedp.Mtb}}, \code{\link{closedp.mX}}, \code{\link{closedp.h}}}
\examples{
data(hare)
hare.closedp<-closedp(hare)
hare.closedp
boxplot(hare.closedp)

data(mvole)
period3<-mvole[,11:15]
closedp(period3)
}
\keyword{models}
