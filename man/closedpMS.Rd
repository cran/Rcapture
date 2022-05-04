\name{closedpMS}
\alias{closedpMS}
\alias{closedpMS.t}
\alias{print.closedpMS}
\alias{plot.closedpMS}
\alias{getAllModels}
\title{ Hierarchical Loglinear Model Selection in Closed Population Capture-Recapture Experiments }
\description{
   The \code{closedpMS.t} function fits every possible hierarchical loglinear models for a given
   closed population capture-recapture data set, under the constraints set by the given \code{maxorder} and \code{forced}
   arguments. Parameters for heterogeneity in capture probabilites among units can be added to the models.
   
   The \code{getAllModels} function lists every possible hierarchical loglinear models for a certain number of capture occasions \code{t}, under the constraints set by the given \code{maxorder} and \code{forced} arguments.
}
\usage{
closedpMS.t(X, dfreq = FALSE, h = NULL, h.control = list(), 
            maxorder = t - 1, forced = 1:t, stopiflong = TRUE, \dots)
			
\method{print}{closedpMS}(x, \dots)

\method{plot}{closedpMS}(x, main="Models comparison based on BIC", omitOutliers = TRUE, \dots)

getAllModels(t, maxorder = t - 1, forced = 1:t, stopiflong = TRUE) 
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} for a 
            description of the accepted formats).  }
  \item{dfreq}{ A logical. By default \code{FALSE}, which means that \code{X} has one row per unit. 
                If \code{TRUE}, it indicates that the matrix \code{X} contains frequencies in its last column.}
  \item{h}{ A character string ("LB", "Chao", "Poisson", "Darroch", "Gamma" or "Normal") or a 
            numerical \code{R} function specifying the form of the column(s) for heterogeneity 
            in the design matrix.\cr
            "LB" and "Chao" represents Chao's lower bound model (the default),\cr 
            "Poisson" represents the function \eqn{f(k)=theta^k-1},\cr 
            "Darroch" represents the function \eqn{f(k)=k^2/2},\cr 
            and "Gamma" represents the function \eqn{f(k)=-\log(theta + k) + \log(theta)}{f(k)=-log(theta + k) + log(theta)}\cr 
            where \eqn{k} is the number of captures and \eqn{theta} is a parameter specified with 
            the \code{h.control} element \code{theta}.\cr
            "Normal" represents heterogeneity 
            modeled with a normal mixing distribution, as presented in Rivest (2011).\cr
            If an \code{R} function is given, it must be the implementation of any convex 
            mathematical function \eqn{f(k)} (it must have only one argument).}
  \item{h.control }{ A list of elements to control the heterogeneous part of the model, if any.\cr
  For a Poisson or Gamma heterogeneous model:
\describe{
  \item{\code{theta}}{ The value of the parameter \eqn{theta} in \eqn{f(k)=theta^k-1} for the Poisson model
                (the default value is 2) and in \eqn{f(k)=-\log(theta + k) + \log(theta)}{f(k)=-log(theta + k) + log(theta)} 
                for the Gamma model (the default value is 3.5). }
}
For a Chao's lower bound heterogeneous model:
\describe{
  \item{\code{neg}}{ If this option is set to \code{TRUE} (the default), negative eta parameters in 
              Chao's lower bound models are set to zero. }
}
For a Normal heterogeneous model:
\describe{       
  \item{\code{initsig}}{ Initial value for the sigma parameter of the normal mixing distribution. 
                  The default is \code{0.2}. }
  \item{\code{method}}{ The method to be used by \code{\link{optim}}. The default is \code{"BFGS"}.\cr }
}
}
  \item{maxorder}{ A numeric specifying the higher order accepted for the terms in the models.
                   It can take a value between 1 and \code{t - 1} (the default).}
  \item{forced}{ A vector of the terms forced in the model (by default every 
                 first order terms). \code{NULL} if no terms are forced in the model. }
  \item{stopiflong}{ A logical indicating whether the function execution should be stopped
                     when the number of possibles models is larger than 10 000, that is when
                     the function might be long to run (defaut \code{TRUE}). }
  \item{\dots}{ Further arguments to be passed to \code{\link{glm}}, \code{\link{optim}} 
                or \code{\link{print.default}}.}
  \item{x}{ An object, produced by the \code{closedpMS.t} function, to print or to plot.}
  \item{main}{ A main title for the plot.}
  \item{omitOutliers}{ A logical. If TRUE (the default), models with an outlier abundance or BIC value are removed from the plot. A value is considered an outlier if it is smaller than the first quartile minus 1.5 times the interquartile range, or larger than the third quartile plus 1.5 times the interquartile range.}
  \item{t}{ A numeric specifying the number of capture occasions in the experiments.
            It is deduced from the data set in the \code{closedpMS.t} function.
            It must be an integer between 1 and 9 inclusively.}
}
\details{
HIERARCHICAL LOGLINEAR MODEL NAME SYNTAX

First, a model's term is written using numbers between 1 and 9 to represent the capture occasions it includes (ex.: \code{134} represents the three-way interaction \code{c1:c3:c4}). This syntax limits the maximal number of capture occation to 9. This is not a problem since from 6 capture occasions upwards, the number of hierarchical models becomes very large and difficult to manage.

A hierarchical model name is a list of the model's terms at the top of the hierarchies in the model. These terms are separated by commas, without spaces. They are surronded by brackets. For example, \code{"[123,34,5]"} is the name of the model\cr
\code{~ 1 + c1 + c2 + c3 + c4 + c5 + c1:c2 + c1:c3 + c2:c3 + c3:c4 + c1:c2:c3}.
}
\value{
\code{getAllModels} returns a caracter vector with the models names.\cr

\code{closedpMS.t} returns a list
with the following elements:
\item{n}{ The number of captured units}
\item{t}{ The total number of capture occasions in the data matrix \code{X}.
          When captures occur in continuous time (input argument \code{t=Inf}),
          this output value is the larger number of captures observed.}
\item{results}{ A table containing, for every fitted model:
    \describe{
    \item{\code{abundance}}{  the estimated population size,}
    \item{\code{stderr}}{  the standard error of the estimated population size,}
    \item{\code{deviance}}{  the model's deviance,}
    \item{\code{df}}{  the number of degrees of freedom,}
    \item{\code{AIC}}{  the Akaike's information criterion,}
    \item{\code{BIC}}{  the bayesian information criterion,}
    \item{\code{infoFit}}{  a numerical code giving information about error or warnings 
          encountered when fitting the model (see \code{\link{Rcapture-package}} for details).}
    }}
\item{bias}{ A vector, the asymptotic bias of the estimated population size for every fitted model.}
\item{fit.err}{ A list of character string vectors. If an error occurs
                while fitting a model (with \code{\link{glm}} or \code{\link{optim}}), a copy of the 
                error message is stored in \code{fit.err$mname} where \code{mname} is the 
                name of the hierarchical loglinear model. A \code{NULL} list element means  
                that no error occured for the considered model. }            
\item{fit.warn}{ A list of character string vectors. If warnings are generated
                 while fitting a model (with \code{\link{glm}} or \code{\link{optim}}), a copy of these 
                 warnings are stored in \code{fit.warn$mname} where \code{mname} is a the 
                 name of the hierarchical model. A \code{NULL} list element means  
                 that no warnings were generated for the considered model. }
\item{neg.eta}{ For Chao's lower bound models only: the position of the eta parameters set 
                to zero in the loglinear parameter vector, if any. } 
}
\note{ For normal heterogeneous models, \code{closedpMS.t} uses \code{\link{optim}}
       from the \pkg{stats} package. Otherwise, models are fitted with \code{\link{glm}},
       also from the \pkg{stats} package. }
\references{ 
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \doi{10.18637/jss.v019.i05}.

Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics}, \bold{43}(4), 783--791.

Cormack, R. M. (1992) Interval estimation for mark-recapture studies of closed populations. \emph{Biometrics}, \bold{48}, 567--576.

Rivest, L.P. (2011) A lower bound model for multiple record systems estimation with heterogeneous catchability. 
\emph{The International Journal of Biostatistics}, \bold{7}(1), Article 23.

Rivest, L.P. and Baillargeon, S. (2007) Applications and extensions of Chao's moment estimator for the size of a closed population. \emph{Biometrics}, \bold{63}(4), 999--1006.
}
\author{ Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and Sophie Baillargeon }

\seealso{\code{\link{closedp}},  \code{\link{closedp.Mtb}}}
\examples{
# The lesbian data set contains 4 capture occasions.
# By default, closedpMB.t fits the 113 following models:

getAllModels(4)

closedpMS.t(lesbian, dfreq = TRUE)

# We could reduce the number of models by omitting
# those with triple interactions.

closedpMS.t(lesbian, dfreq = TRUE, maxorder = 2)

# Models with heterogeneity fits better.

Darr <- closedpMS.t(lesbian, dfreq = TRUE, h = "Darroch")
Darr

# The plot method allows the visualization of the results
# from models fitted by closedpMS.t().

plot(Darr)

# According to the BIC, the best heterogeneous Darroch model
# for this data set contains the double interactions 12, 13, 14.
# Here is the profile likelihood confidence interval for the
# abundance estimation from this model.

closedpCI.t(lesbian, dfreq = TRUE, mX = "[12,13,14]", h = "Darroch")
}
\keyword{ models }
