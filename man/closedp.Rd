\name{closedp}
\alias{closedp}
\alias{closedp.t}
\alias{closedp.0}
\alias{print.closedp}
\alias{boxplot.closedp}
\alias{plot.closedp}
\title{ Loglinear Models for Closed Population Capture-Recapture Experiments }
\description{
The functions \code{closedp.t} and \code{closedp.0} fit various loglinear models for closed populations in 
capture-recapture experiments. For back compatibility, \code{closedp.t} is also named \code{closedp}.

\code{closedp.t} fits more models than \code{closedp.0} but for data set with more than 20 capture occasions, the function might fail. However, \code{closedp.0} works with fairly large data sets (see \bold{Details}).     
}
\usage{
closedp(X, dfreq=FALSE, neg=TRUE, \dots)
closedp.t(X, dfreq=FALSE, neg=TRUE, \dots)

closedp.0(X, dfreq=FALSE, dtype=c("hist","nbcap"), t=NULL, t0=NULL, 
          neg=TRUE, \dots)

\method{print}{closedp}(x, \dots)

\method{boxplot}{closedp}(x, main="Boxplots of Pearson Residuals", \dots)

\method{plot}{closedp}(x, main="Residual plots for some heterogeneity models", \dots)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} 
            for a description of the accepted formats).  }
  \item{dfreq}{ A logical. By default \code{FALSE}, which means that \code{X} has one row per unit. 
                If \code{TRUE}, it indicates that the matrix \code{X} contains frequencies in its last column.}
  \item{dtype}{ A characters string, either \code{"hist"} or \code{"nbcap"}, to specify the type of data. 
                \code{"hist"}, the default, means that \code{X} contains complete observed capture histories. 
                \code{"nbcap"} means that \code{X} contains numbers of captures (see \code{\link{Rcapture-package}} 
                for details on data formats). }
  \item{t}{ Requested only if \code{dtype="nbcap"}. A numeric specifying the total number of 
            capture occasions in the experiment. For \code{closedp.0}, the value \code{t=Inf}
            is accepted. It indicates that captures occur in continuous time 
            (see \code{\link{Rcapture-package}}).}
  \item{t0}{ A numeric. Models are fitted considering only the frequencies of units captured 1 to \code{t0} times. 
             By default, if \code{t} is not equal to \code{Inf}, \code{t0=t}. When \code{t=Inf}, the default value
             of \code{t0} is the larger number of captures observed. }
  \item{neg}{ If this option is set to TRUE, negative eta parameters in Chao's lower 
              bound models are set to zero (see \bold{Details}). }
  \item{\dots}{ Further arguments to be passed to \code{\link{glm}}, \code{\link{print.default}},
                \code{\link{boxplot.default}} or \code{\link{plot.default}}.}
  \item{x}{ An object, produced by a \code{closedp} function, to print or to plot.}
  \item{main}{ A main title for the plot.}
  }
\details{
   \code{closedp.t} fits models M0, Mt, Mh Chao (LB), Mh Poisson2, Mh Darroch, Mh Gamma3.5, 
   Mth Chao (LB), Mth Poisson2, Mth Darroch, Mth Gamma3.5, Mb and Mbh. \code{closedp.0} fits 
   only models M0, Mh Chao (LB), Mh Poisson2, Mh Darroch and Mh Gamma3.5. However, 
   \code{closedp.0} can be used with larger data sets than \code{closedp.t}. 
   This is explained by the fact that \code{closedp.t} fits models using the frequencies of 
   the observable capture histories (vector of size \eqn{2^t-1}), whereas \code{closedp.0} 
   uses the numbers of units captured i times, for \eqn{i=1,\ldots,t} (vector of size \eqn{t}). 
   See \code{\link{Rcapture-package}} for more details about the distinction between \code{.t}
   and \code{.0} functions.
    
  Multinomial profile confidence intervals for the abundance are constructed by \code{\link{closedpCI.t}} 
  and \code{\link{closedpCI.0}}.

  To calculate bias corrected abundance estimates, use the \code{\link{closedp.bc}} function.\cr

  CHAO'S LOWER BOUND MODEL
        
  Chao's (or LB) models estimate a lower bound for the abundance, both with a time 
  effect (Mth Chao) and without one (Mh Chao). The estimate obtained under Mh Chao is Chao's (1987) 
  moment estimator. Rivest and Baillargeon (2007) exhibit a loglinear model underlying this 
  estimator and provide a generalization to Mth. For these two models, a small deviance means 
  that there is an heterogeneity in capture probabilities; it does not mean that the lower 
  bound estimate is unbiased. To test whether a certain model for heterogeneity is adequate,
  one can conduct a likelihood ratio test by subtracting the deviance of Chao's model
  to the deviance of the heterogeneous model under study. If this heterogeneous model
  includes a time effect, it must be compared to model Mth Chao. If it does not include a time
  effect, it must be compared to model Mh Chao. 
  Under the null hypothesis of equivalence between the two models, the difference of deviances 
  follows a chi-square distribution with degrees of freedom equal to the difference between 
  the models' degrees of freedom.
       
  Chao's lower bound models contain \eqn{t-2} parameters, called
  eta parameters, for the heterogeneity. These parameters should theoretically be greater 
  or equal to zero (see Rivest and Baillargeon (2007)). When the argument \code{neg} is set 
  to \code{TRUE} (the default), negative eta parameters are set to zero (to do so, columns are
  removed from the design matrix of the model). Degrees of freedom of Chao's model increase 
  when eta parameters are set to zero.\cr

OTHER MODELS FOR HETEROGENEITY

  Other models for heterogeneity are defined as follows :
  \tabular{ll}{
  Model \tab Column for heterogeneity in the design matrix \cr
  \tab \cr
  Poisson2 \tab \eqn{2^k-1} \cr
  Darroch  \tab \eqn{k^2/2} \cr
  Gamma3.5 \tab \eqn{-\log(3.5 + k) + \log(3.5)}{-log(3.5 + k) + log(3.5)}
  }
  where \eqn{k} is the number of captures. Poisson and Gamma models with alternative to the 
  parameter defaults values 2 and 3.5 can be fitted with the \code{\link{closedpCI.t}} and 
  \code{\link{closedpCI.0}} functions.

  Darroch's models for Mh and Mth are considered by Darroch et al. (1993) and Agresti (1994). 
  Poisson and Gamma models are discussed in Rivest and Baillargeon (2007). Poisson models 
  typically yield smaller corrections for heterogeneity than Darroch's model since the capture 
  probabilities are bounded from below under these models. On the other hand, Gamma models 
  can lead to very large estimators of abundance. We suggest considering this estimator only in 
  experiments where very small capture probabilities are likely.\cr

  
PLOT METHODS AND FUNCTIONS
  
  The \code{boxplot.closedp} function produces boxplots of the Pearson residuals of the fitted loglinear models that converged.

  The \code{plot.closedp} function produces scatterplots of the Pearson residuals in terms of \eqn{f_i}{fi} 
  (number of units captured i times) for the heterogeneous models Mh Poisson2, Mh Darroch and Mh Gamma3.5 if they converged.  
}
\value{
\item{n}{ The number of captured units.}
\item{t}{ The total number of capture occasions in the data matrix \code{X}.
          When captures occur in continuous time (input argument \code{t=Inf}),
          this output value is the larger number of captures observed.}
\item{t0}{ For \code{closedp.0} only: the value of the argument \code{t0} used in the computations.}
\item{results}{ A table containing, for every fitted model:
    \describe{
    \item{\code{abundance}}{:  the estimated population size,}
    \item{\code{stderr}}{:  the standard error of the estimated population size,}
    \item{\code{deviance}}{:  the model's deviance,}
    \item{\code{df}}{:  the number of degrees of freedom,}
    \item{\code{AIC}}{:  the Akaike's information criterion,}
    \item{\code{BIC}}{:  the bayesian information criterion,}
    \item{\code{infoFit}}{:  a numerical code giving information about error or warnings 
          encountered when fitting the model (see \code{\link{Rcapture-package}} for details).}
    }}
\item{bias}{ A vector, the asymptotic bias of the estimated population size for every fitted model.}
\item{glm}{ A list of the 'glm' objects obtained from fitting models.} 
\item{glm.err}{ A list of character string vectors. If the \code{\link{glm}} function generates
                an error when fitting a model, a copy of the error message is
                 stored in \code{glm.err$mname} where \code{mname} is a short name to identify
                 the model. A \code{NULL} list element means that no error occured for the considered model.} 
\item{glm.warn}{ A list of character string vectors. If the \code{\link{glm}} function generates
                 one or more warnings when fitting a model, a copy of these warnings are
                 stored in \code{glm.warn$mname} where \code{mname} is a short name to identify
                 the model. A \code{NULL} list element means that \code{\link{glm}} did not 
                 produce any warnings for the considered model. }
\item{parameters}{ Capture-recapture parameters estimates. It contains N, the estimated population size, 
                   and p or \eqn{p_1}{p1} to \eqn{p_t}{pt} defined as follows for the different models :
    \describe{
    \item{}{M0: the capture probability at any capture occasion,}
    \item{}{Mt: the capture probabilities for each capture occasion,}
    \item{}{Mh models: the average probability of capture,} 
    \item{}{Mth models: the average probabilities of capture for each occasion,}
    \item{}{Mb and Mbh: the probability of first capture at any capture occasion.} 
    }
    For models Mb and Mbh, it also contains c, the recapture probability at any capture occasion.}
\item{neg.eta}{ The position of the eta parameters set to zero in the loglinear parameter 
                vector of models MhC and MthC. \code{NULL} if no parameters are set to zero.}
\item{X}{ A copy of the data given as input in the function call.}
\item{dfreq}{ A copy of the \code{dfreq} argument given in the function call.}
}
\note{ This function uses the \code{\link{glm}} function of the \pkg{stats} package. }
\references{ 
Agresti, A. (1994) Simple capture-recapture models permitting unequal catchability and variable sampling effort. \emph{Biometrics}, \bold{50}, 494--500.

Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.

Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchabililty.  \emph{ Biometrics}, \bold{45}, 427--438.

Darroch, S.E.,  Fienberg, G.,  Glonek, B. and Junker, B. (1993) A three sample multiple capture-recapture approach to the census population estimation with heterogeneous catchability. \emph{Journal of the American Statistical Association}, \bold{88}, 1137--1148.

Rivest, L.P. and Levesque, T. (2001) Improved loglinear model estimators of abundance in capture-recapture experiments. \emph{Canadian Journal of Statistics}, 29, 555--572.

Rivest, L.P. and Baillargeon, S. (2007) Applications and extensions of Chao's moment estimator for the size of a closed population. \emph{Biometrics}, \bold{63}(4), 999--1006.

Seber, G.A.F. (1982) The Estimation of Animal Abundance and Related Parameters, 2nd edition, New York: Macmillan.
}
\author{ Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and Sophie Baillargeon }

\seealso{\code{\link{closedpCI.t}}, \code{\link{closedpCI.0}}, \code{\link{closedp.bc}}, \code{\link{closedp.Mtb}}, \code{\link{closedpMS.t}}, \code{\link{uifit}}. }
\examples{
data(hare)
hare.closedp <- closedp.t(hare)
hare.closedp
boxplot(hare.closedp)

data(mvole)
period3 <- mvole[, 11:15]
closedp.t(period3)

data(BBS2001)
BBS.closedp <- closedp.0(BBS2001, dfreq = TRUE, dtype = "nbcap", t = 50, t0 = 20)
BBS.closedp
plot(BBS.closedp)

### Seber (1982) p.107
# When there is 2 capture occasions, the heterogeneity models cannot be fitted
X <- matrix(c(1,1,167,1,0,781,0,1,254), byrow = TRUE, ncol = 3)
closedp.t(X, dfreq = TRUE)

### Example of captures in continuous time
### Illegal immigrants data
data(ill)
closedp.0(ill, dtype = "nbcap", dfreq = TRUE, t = Inf)
}
\keyword{models}
