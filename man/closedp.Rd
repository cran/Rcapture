\name{closedp}
\alias{closedp}
\alias{closedp.t}
\alias{closedp.0}
\alias{print.closedp}
\alias{boxplot.closedp}
\alias{plot.closedp}
\title{ Loglinear Models for Closed Population Capture-Recapture Experiments }
\description{
The functions \code{closedp.t} and \code{closedp.0} fit various loglinear models for closed populations in capture-recapture experiments. For back compatibility, \code{closedp.t} is also named \code{closedp}.

\code{closedp.t} fits more models than \code{closedp.0} but for data set with more than 20 capture occasions, the function migth fail. However, \code{closedp.0} works with fairly large data sets (see Details).     
}
\usage{
closedp(X, dfreq=FALSE, neg=TRUE, trace=FALSE)
closedp.t(X, dfreq=FALSE, neg=TRUE, trace=FALSE)

closedp.0(X, dfreq=FALSE, dtype=c("hist","nbcap"), t, t0=t, 
          neg=TRUE, trace=FALSE)

\method{print}{closedp}(x, \dots)

\method{boxplot}{closedp}(x, main="Boxplots of Pearson Residuals", \dots)

\method{plot}{closedp}(x, main="Residual plots for some heterogeneity models", \dots)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} for a description of the accepted formats).  }
  \item{dfreq}{ A logical. By default FALSE, which means that \code{X} has one row per unit. If TRUE, it indicates that the matrix \code{X} contains frequencies in its last column.}
  \item{dtype}{ A characters string, either "hist" or "nbcap", to specify the type of data. "hist", the default, means that \code{X} contains complete observed capture histories. "nbcap" means that \code{X} contains numbers of captures (see \code{\link{Rcapture-package}} for details on data formats). }
  \item{t}{ Requested only if \code{dtype="nbcap"}. A numeric specifying the total number of capture occasions in the experiment.}
  \item{t0}{ A numeric. Models are fitted considering only the frequencies of units captured 1 to \code{t0} times. By default \code{t0=t}. }
  \item{neg}{ If this option is set to TRUE, negative eta parameters in Chao's models are set to zero. }
  \item{trace}{ A logical, by default FALSE. If set to TRUE, a note is printed for each model while the function runs. It is useful to identify which model is associated to a warning. }
  \item{x}{ An object, produced by a \code{closedp} function, to print or to plot.}
  \item{main}{ A main title for the plot}
  \item{\dots}{ Further arguments to be passed to methods (see \code{print.default}, \code{boxplot.default} and \code{plot.default}).}
  }
\details{
   \code{closedp.t} fits models M0, Mt, Mh Chao, Mh Poisson2, Mh Darroch, Mh Gamma3.5, Mth Chao, Mth Poisson2, Mth Darroch, Mth Gamma3.5, Mb and Mbh. \code{closedp.0} fits only models M0, Mh Chao, Mh Poisson2, Mh Darroch and Mh Gamma3.5. However, \code{closedp.0} can be used with larger data sets than \code{closedp.t}. This is explained by the fact that \code{closedp.t} fits models using the frequencies of the observable capture histories (vector of size \eqn{2^t-1}), whereas \code{closedp.0} uses the numbers of units captured i times, for \eqn{i=1,\ldots,t} (vector of size \eqn{t}).
   
   \code{closedp.0} has an additional argument \code{t0} which gives to the numbers of units caught more than \code{t0} times their own parameters in the loglinear model. For example, the model for Mh Gamma3.5 has \eqn{3+t-t_0}{3+t-t0} parameters. This means that \code{closedp.0} fits models considering only the frequencies of units captured 1 to \code{t0} times.     
   
  Chao's models estimate a lower bound for the abundance, both with a time effect (Mth Chao) and without (Mh Chao). The estimate obtained under Mh Chao is Chao's (1987) moment estimator. Rivest and Baillargeon (2007) exhibit a loglinear model underlying this estimator and provide a generalization to Mth. For these two models, a small deviance means that there is an heterogeneity in capture probabilities; it does not mean that the lower bound estimates are unbiased.

  Other models for heterogeneity are defined as follows :
  \tabular{ll}{
  Model \tab Colomn for heterogeneity in the design matrix \cr
  \tab \cr
  Poisson2 \tab \eqn{2^k-1} \cr
  Darroch  \tab \eqn{k^2/2} \cr
  Gamma3.5 \tab \eqn{-log(3.5 + k) + log(3.5)}
  }
  where \eqn{k} is the number of captures. Poisson and Gamma models with parameter values different than 2 and 3.5 respectively can be fitted with the \code{\link{closedpCI.t}} and \code{\link{closedpCI.0}} functions.

  Darroch's models for Mh and Mth are considered by Darroch et al. (1993) and Agresti (1994). Poisson and Gamma models are discussed in Rivest and Baillargeon (2007). Poisson models typically yield smaller corrections for heterogeneity than Darroch's model since the capture probabilities are bounded from below under these models. On the other hand, Gamma models can lead to very large estimators of abundance. We suggest considering this estimator in experiments where very small capture probabilities are likely.

  When the variance of an abundance estimate is large, it is useful to use the \code{\link{closedpCI.t}} or \code{\link{closedpCI.0}} function to construct a profile condifence interval for this abundance.
  
  The \code{boxplot.closedp} function produces boxplots of the Pearson residuals of the fitted loglinear models that converged.

  The \code{plot.closedp} function produces scatterplots of the Pearson residuals in terms of \eqn{f_i}{fi} (number of units captured i times) for the heterogeneous models Mh Poisson2, Mh Darroch and Mh Gamma3.5 if they converged.  
  
  To calculate bias corrected abundance estimates, use the \code{\link{closedp.bc}} function.
}
\value{
\item{n}{ The number of captured units.}
\item{t}{ The total number of capture occasions in the experiment.}
\item{results}{ A table containing, for every fitted model, the estimated population size, the standard error of estimation, the deviance, the number of degrees of freedom and the Akaike criteria.}
\item{converge}{ A logical vector indicating whether or not the fitted models converged.}
\item{glm}{ A list of the 'glm' objects obtained from fitting models.} 
\item{parameters}{ Capture-recapture parameters estimates. It contains N, the estimated population size, and p or \eqn{p_1}{p1} to \eqn{p_t}{pt} defined as follows for the different models :
\tabular{ll}{
M0 \tab the capture probability at any capture occasion \cr
Mt \tab the capture probabilities for each capture occasion \cr
Mh models \tab the average probability of capture \cr 
Mth models \tab the average probabilities of capture for each occasion \cr
Mb and Mbh \tab the probability of first capture at any capture occasion 
}
For models Mb and Mbh, it also contains c, the recapture probability at any capture occation.}
\item{neg.eta}{ The position of the eta parameters set to zero in the loglinear parameter vector of model MhC and MthC.}
\item{X}{ A copy of the data given as input in the function call.}
\item{dfreq}{ A copy of the \code{dfreq} argument given in the function call.}
\item{t0}{ For \code{closedp.0} only, a copy of the \code{t0} argument given in the function call.}
}
\references{ 
Agresti, A. (1994) Simple capture-recapture models permitting unequal catchability and variable sampling effort. \emph{Biometrics}, \bold{50}, 494--500.

Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.

Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchabililty.  \emph{ Biometrics}, \bold{45}, 427--438.

Darroch, S.E.,  Fienberg, G.,  Glonek, B. and Junker, B. (1993) A three sample multiple capture-recapture approach to the census population estimation with heterogeneous catchability. \emph{Journal of the American Statistical Association}, \bold{88}, 1137--1148.

Rivest, L.P. and Levesque, T. (2001) Improved log-linear model estimators of abundance in capture-recapture experiments. \emph{Canadian Journal of Statistics}, 29, 555--572.

Rivest, L.P. and Baillargeon, S. (2007) Applications and extensions of Chao's moment estimator for the size of a closed population. \emph{Biometrics}, \bold{63}(4), 999--1006.

Seber, G.A.F. (1982) The Estimation of Animal Abundance and Related Parameters, 2nd edition, New York: Macmillan.
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\note{ This function uses the \code{glm} function of the \code{stats} package. }
\seealso{\code{\link{closedpCI.t}}, \code{\link{closedpCI.0}}, \code{\link{closedp.bc}}, \code{\link{closedp.Mtb}}, \code{\link{uifit}}. }
\examples{
data(hare)
hare.closedp<-closedp.t(hare)
hare.closedp
boxplot(hare.closedp)

data(mvole)
period3<-mvole[,11:15]
closedp.t(period3)

data(BBS2001)
BBS.closedp<-closedp.0(BBS2001,dfreq=TRUE,dtype="nbcap",t=50,t0=20)
BBS.closedp
plot(BBS.closedp)

# Seber (1982) p.107
# When there is 2 capture occasions, the heterogeneity models cannot be fitted
X <- matrix(c(1,1,167,1,0,781,0,1,254),byrow=TRUE,ncol=3)
closedp.t(X,dfreq=TRUE)
}
\keyword{models}
