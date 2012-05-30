\name{closedpCI}
\alias{closedpCI}
\alias{closedpCI.t}
\alias{closedpCI.0}
\alias{print.closedpCI}
\alias{plotCI}
\alias{boxplot.closedpCI}
\alias{plot.closedpCI}
\title{ Customization of a Loglinear Model and Confidence Interval for Abundance 
        Estimation in Closed Population Capture-Recapture Experiments }
\description{
   The \code{closedpCI.t} and \code{closedpCI.0} functions fit a loglinear model specified by the user 
   and compute a confidence interval for the abundance estimation. For a normal heterogeneous model,
   a log-transformed confidence interval (Chao 1987) is produced.
   For any other model, the multinomial profile likelihood confidence interval (Cormack 1992) is produced.\cr
   
   The model is identified with the argument \code{m} or \code{mX}.
   For heterogeneous models, the form of the heterogeneity is specified with the arguments 
   \code{h} and \code{h.control}. If \code{h} is given with \code{mX}, heterogeneity is added in \code{mX}.\cr
   
   These functions extend \code{closedp.t} and \code{closedp.0} as 
   they broaden the range of models one can fit and they compute confidence intervals. 
   Unlike the \code{closedp} functions, it fits only one model at a time. 
}
\usage{
closedpCI.t(X, dfreq=FALSE, m=c("M0","Mt","Mh","Mth"), mX=NULL,
            h=NULL, h.control=list(), mname=NULL, alpha=0.05, 
            fmaxSupCL=3, \dots)

closedpCI.0(X, dfreq=FALSE, dtype=c("hist","nbcap"), t=NULL, t0=NULL, 
            m=c("M0","Mh"), mX=NULL, h=NULL, h.control=list(), 
            mname=NULL, alpha=0.05, fmaxSupCL=3, \dots)
			
\method{print}{closedpCI}(x, \dots)

\method{boxplot}{closedpCI}(x, main="Boxplots of Pearson Residuals", \dots)

\method{plot}{closedpCI}(x, main="Scatterplot of Pearson Residuals", \dots)

plotCI(x.closedpCI, main="Profile Likelihood Confidence Interval", \dots)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} for a 
            description of the accepted formats).  }
  \item{dfreq}{ A logical. By default \code{FALSE}, which means that \code{X} has one row per unit. 
                If \code{TRUE}, it indicates that the matrix \code{X} contains frequencies in its last column.}
  \item{dtype}{ A characters string, either \code{"hist"} or \code{"nbcap"}, to specify the type of data. 
                \code{"hist"}, the default, means that \code{X} contains complete observed capture histories. 
                \code{"nbcap"} means that \code{X} contains numbers of captures (see \code{\link{Rcapture-package}} 
                for details on data formats). }
  \item{t}{ Requested only if \code{dtype="nbcap"}. A numeric specifying the total number of 
            capture occasions in the experiment. For \code{closedpCI.0}, the value \code{t=Inf}
            is accepted. It indicates that captures occur in continuous time 
            (see \code{\link{Rcapture-package}}).}
  \item{t0}{ A numeric. Models are fitted considering only the frequencies of units captured 1 to \code{t0} times. 
             By default, if \code{t} is not equal to \code{Inf}, \code{t0=t}. When \code{t=Inf}, the default value
             of \code{t0} is the larger number of captures observed. }
  \item{m}{ A character string indicating the model to fit. For \code{closedpCI.0} it can be either "M0"=M0 model 
            or "Mh"=Mh model. For \code{closedpCI.t} it can also be "Mt"=Mt model or "Mth"=Mth model.}
  \item{mX}{ The design matrix of the loglinear model. By default, the design matrix is built based
             on the \code{m} argument. If a \code{mX} argument is given, it must be a matrix  (or 
             an object that can be coerced to a matrix by \code{as.matrix}). The order of the capture 
             histories in the \code{mX} matrix must be as defined in the \code{histpos.t} or 
             \code{histpos.0} function. For the \code{closedpCI.t} function only, \code{mX}
             can also be an object of class "\code{\link{formula}}". The details of a formula 
             specification are given under \bold{Details}. An intercept is always added to the 
             model. Therefore, \code{mX} must not  
             contain a column of ones if it is a matrix. Also, if the argument \code{h} 
             is not \code{NULL}, one or more columns for heterogeneity are added to the design 
             matrix.}
  \item{h}{ A character string ("Chao", "LB", "Poisson", "Darroch", "Gamma" or "Normal") or a 
            numerical \code{R} function specifying the form of the column(s) for heterogeneity 
            in the design matrix. "Chao" and "LB" represents Chao's lower bound model, 
            "Poisson" represents the function \eqn{f(k)=theta^k-1}, 
            "Darroch" represents the function \eqn{f(k)=k^2/2}, 
            and "Gamma" represents the function \eqn{f(k)=-\log(theta + k) + \log(theta)}{f(k)=-log(theta + k) + log(theta)} 
            where \eqn{k} is the number of captures and \eqn{theta} is a parameter specified with 
            the \code{h.control} element \code{theta}. "Normal" represents heterogeneity 
            modeled with a normal mixing distribution, as presented in Rivest (2011).
            If an \code{R} function is given, it must be the implementation of any convex 
            mathematical function \eqn{f(k)} (it must have only one argument).}
  \item{h.control }{ A list of elements to control the heterogeneous part of the model, if any (see \bold{Details}). }
  \item{mname}{ A character string specifying the name of the customized model. By default,
                it is derived from the arguments specifying the model.}
  \item{alpha}{ A confidence interval with confidence level 1-\code{alpha} is constructed. The value of alpha 
                must be between 0 and 1; the default is 0.05.  }
  \item{fmaxSupCL}{ A numeric. The upper end point of the interval to be searched by \code{uniroot} 
                   to find the upper bound of the multinomial profile likelihood confidence interval 
                   (Cormack 1992) is defined by \code{fmaxSupCL}\eqn{\times \hat{N}}{x \hat{N}}. In 
                   this product, \code{fmaxSupCL} is a factor multiplying \eqn{\hat{N}}, the population 
                   size estimated by Poisson regression. If the upper bound obtained is equal
                   the upper end point \code{fmaxSupCL}\eqn{\times \hat{N}}{x \hat{N}}, the 
                   element \code{SupCL} of the output value \code{CI} will begin with the symbol 
                   \eqn{>}. In that case, one could try to increase the value of \code{fmaxSupCL} 
                   to find the upper bound of the multinomial profile likelihood confidence interval 
                   instead of a minimum value for this upper bound. }
  \item{x}{ An object, produced by a \code{closedpCI} function, to print or to plot.}
  \item{x.closedpCI}{ An object, produced by a \code{closedpCI} function, to produce a plot
                      of the multinomial profile likelihood for \eqn{N}.}
  \item{main}{ A main title for the plot}
  \item{\dots}{ Further arguments to be passed to \code{optim}, \code{print.default}, \code{plot.default} or 
                \code{boxplot.default}.}
  }
\details{ 
The \code{closedpCI.t} function fits models using the frequencies of the observable capture histories 
(vector of size \eqn{2^t-1}), whereas \code{closedpCI.0} uses the number of units capture i times, for 
\eqn{i=1,\ldots,t}{i=1,\ldots,t} (vector of size \eqn{t}). Thus, \code{closedpCI.0} can be used with data 
sets larger than those for \code{closedpCI.t}, but it cannot fit models with a temporal effect. 

These functions do not work for closed population models featuring a behavioral effect, such as Mb and Mbh. 
The abundance estimation is calculated as the number of captured units plus the exponential of the Poisson 
regression intercept. However, models with a behavioral effect can by fitted with \code{\link{closedp.t}} 
(Mb and Mbh), \code{\link{closedp.Mtb}} and \code{\link{closedp.bc}}.\cr

CHAO'S LOWER BOUND MODELS

  Chao's lower bound models estimate a lower bound for the abundance. Rivest (2011) 
  presents a generalized loglinear model underlying this estimator. To test whether a 
  certain model for heterogeneity is adequate,
  one can conduct a likelihood ratio test by subtracting the deviance of a Chao's lower 
  bound model to the deviance of the heterogeneous model under study. The two models should
  have the same \code{mX} argument.  
  Under the null hypothesis of equivalence between the two models, the difference of deviances 
  follows a chi-square distribution with degrees of freedom equal to the difference between 
  the models' degrees of freedom.
       
  A Chao's lower bound model contains \eqn{t-2} parameters, called
  eta parameters, for the heterogeneity. These parameters should theoretically be greater 
  or equal to zero (see Rivest and Baillargeon (2007)). When the element \code{neg}
  of the argument \code{h.control} is set 
  to TRUE (the default), negative eta parameters are set to zero (to do so, columns are
  removed from the design matrix of the model). Consequently, heterogeneous models are no 
  longer particular cases of 
  Chao's lower bound model. Therefore, likelihood ratio tests cannot be conducted 
  anymore to test whether a chosen heterogeneous model is adequate. Also, the degrees 
  of freedom of Chao's model increase when eta parameters are set to zero.\cr

ARGUMENT \code{mX} : formula specification\cr
For the \code{closedpCI.t} function, \code{mX}
can be an object of class "\code{\link{formula}}". The only accepted variables 
in this formula are \code{c1} to \code{ct}. The variable \code{ci} represents 
a capture indicator (1 for a capture, 0 otherwise) for the \eqn{i}th capture occasions. 
Also, the formula must not contain a response variable since
it is only used to construct the design matrix of the model. 
For example, if \code{t=3}, the Mt model is fitted if 
\code{mX = ~ .} or \code{mX = ~ c1 + c2 + c3}. The symbol \code{.} in this formula
is a shortcut for \code{c1} +  \code{c2} + ... + \code{ct}. Formula \code{mX} arguments
facilitate the addition of interactions between capture occasions in the model. For
example, if \code{t=3}, the Mt model with an interaction between the first and
the second capture occasion is fitted if \code{mX = ~ . + c1:c2}. 
See \code{\link{formula}} for more details of allowed formulae.\cr

ARGUMENT \code{h.control}\cr
The \code{h.control} argument is a list to supply any of the following elements to control the 
heterogeneous part of the model, if any.

For a Poisson or Gamma heterogeneous model:
\describe{
  \item{theta}{ The value of the parameter \eqn{theta} in \eqn{f(k)=theta^k-1} for the Poisson model
                (the default value is 2) and in \eqn{f(k)=-\log(theta + k) + \log(theta)}{f(k)=-log(theta + k) + log(theta)} 
                for the Gamma model (the default value is 3.5). }
}
For a Chao's lower bound heterogeneous model:
\describe{
  \item{neg}{ If this option is set to \code{TRUE} (the default), negative eta parameters in 
              Chao's lower bound models are set to zero. }
}
For a Normal heterogeneous model:
\describe{       
  \item{initcoef}{ Initial values for the loglinear coefficients given to \code{\link{optim}}. 
                   The default values are taken from Darroch's model. }
  \item{initsig}{ Initial value for the sigma parameter of the normal mixing distribution. 
                  The default is \code{0.2}. }
  \item{method}{ The method to be used by \code{\link{optim}}. The default is \code{"BFGS"}.\cr }
}

PLOT METHODS AND FUNCTIONS
  
  The \code{boxplot.closedpCI} function produces a boxplot of the Pearson residuals of the customized model.

  The \code{plot.closedpCI} function traces the scatterplot of the Pearson residuals in terms of \eqn{f_i}{fi} 
  (number of units captured i times) for the customized model.  
  
  The \code{plotCI} function produces a plot of the multinomial profile likelihood for \eqn{N}. 
  The value of N maximizing the profile likelihood and the bounds of the confidence interval are identified.
  }
\value{
\item{n}{ The number of captured units}
\item{t}{ The number of capture occasions in the data matrix \code{X}}
\item{t0}{ For \code{closedpCI.0} only: the value of the argument \code{t0} used in the computations.}
\item{results}{ A table containing the estimated population size and its standard error,
                the deviance, the number of degrees of freedom and the Akaike's information criterion. If the name
                of the model is followed by ** in this table, it means that the model did not 
                converge. Therefore, the estimated population size for this model is 
                questionable.}
\item{fit}{ The 'glm' object obtained from fitting the model except for
            normal heterogeneous models (\code{h="Normal"}). These models are not fitted 
            with \code{glm.fit}. For them, \code{fit} is a list with the following elements:\cr
            \code{parameters}: The matrix of parameters (loglinear coefficients + sigma parameter) 
                               estimates with their standard errors.\cr
            \code{varcov}: The estimated variance-covariance matrix of the estimated parameters.\cr
            \code{y}: The y vector used to fit the model.\cr
            \code{fitted.values}: The model fitted values.\cr 
            \code{initparam}: The initial values for the parameters (loglinear coefficients + 
                              sigma parameter) used by \code{\link{optim}}.\cr
            \code{optim.out}: The output produced by \code{\link{optim}}.}            
\item{glm.warn}{ Not produced for normal heterogeneous models (\code{h="Normal"}): 
                 A vector of character strings. If the \code{\link{glm.fit}} function generates
                 one or more warnings when fitting the model, a copy of these warnings are
                 stored in \code{glm.warn}. \code{NULL} if \code{\link{glm.fit}} did not produce
                 any warnings. }
\item{neg.eta}{ For Chao's lower bound models only: the position of the eta parameters set 
                to zero in the loglinear parameter vector, if any. If not \code{NULL}, the deviance
                and degrees of freedom of the fitted Chao's lower bound model cannot be
                used to conduct a likelihood ratio test to investigate whether a particular 
                heterogeneous model is adequate (see \bold{Details}). } 
\item{CI}{ Not produced for normal heterogeneous models (\code{h="Normal"}): A table containing 
           the abundance estimation and its multinomial profile likelihood confidence interval.}
\item{alpha}{ 1-the confidence level of the interval.}
\item{NCI}{ Not produced for normal heterogeneous models (\code{h="Normal"}): 
            The x-coordinates for \code{plot.closedpCI.t}.}
\item{loglikCI}{ Not produced for normal heterogeneous models (\code{h="Normal"}): 
                 The y-coordinates for \code{plot.closedpCI.t}.} 
}
\note{ For normal heterogeneous models, the \code{closedpCI} functions use \code{\link{optim}}
       from the \pkg{stats} package. Otherwise, models are fitted with \code{\link{glm.fit}}
       and the code to compute the multinomial profile likelihood confidence interval 
       calls the functions \code{\link{glm.fit}}, \code{\link{optimize}} and \code{\link{uniroot}},
       also from the \pkg{stats} package. }
\references{ 
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.

Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics}, \bold{43}(4), 783--791.

Cormack, R. M. (1992) Interval estimation for mark-recapture studies of closed populations. \emph{Biometrics}, \bold{48}, 567--576.

Rivest, L.P. (2011) A lower bound model for multiple record systems estimation with heterogeneous catchability. 
\emph{The International Journal of Biostatistics}, \bold{7}(1), Article 23.

Rivest, L.P. and Baillargeon, S. (2007) Applications and extensions of Chao's moment estimator for the size of a closed population. \emph{Biometrics}, \bold{63}(4), 999--1006.
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\seealso{\code{\link{closedp}},  \code{\link{closedp.Mtb}}}
\examples{
data(hare)
CI<-closedpCI.t(hare, m = "Mth", h = "Poisson", h.control = list(theta = 2))
CI
plotCI(CI)

data(HIV)
mat<-histpos.t(4)
mX2<-cbind(mat,mat[,1]*mat[,2])
closedpCI.t(HIV, dfreq = TRUE, mX = mX2, mname = "Mt interaction 1,2")
# which can be obtained more conveniently with
closedpCI.t(HIV, dfreq = TRUE, mX = ~ . + c1:c2, mname = "Mt interaction 1,2")

data(BBS2001)
CI0<-closedpCI.0(BBS2001, dfreq = TRUE, dtype = "nbcap", t = 50, t0 = 20,
                 m = "Mh", h = "Gamma", h.control = list(theta = 3.5))
CI0
plot(CI0)
plotCI(CI0)

### As an alternative to a gamma model, one can fit a negative Poisson model.
### It is appropriate in experiments where very small capture probabilities
### are likely. It can lead to very large estimators of abundance. 
data(mvole)
period3 <- mvole[, 11:15]
psi <- function(x) { 0.5^x - 1 }
closedpCI.t(period3, m = "Mh", h = psi)

### Example of normal heterogeneous models
### diabetes data of Bruno et al. (1994)
histpos <- histpos.t(4)
diabetes <- cbind(histpos, c(58,157,18,104,46,650,12,709,14,20,7,74,8,182,10))
# chosen interaction set I in Rivest (2011)
closedpCI.t(X=diabetes, dfreq=TRUE, mX= ~ . + c1:c3 + c2:c4 + c3:c4, 
            h="Normal", mname="Mth normal with I")

### Example of captures in continuous time
### Illegal immigrants data
data(ill)
closedpCI.0(ill, dtype="nbcap", dfreq=TRUE, t=Inf, m="Mh", h="LB")
}
\keyword{ models }
