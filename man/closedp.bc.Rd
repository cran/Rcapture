\name{closedp.bc}
\alias{closedp.bc}
\alias{print.closedp.bc}
\title{ Bias Correction for Closed Population Loglinear Models }
\description{
  This function applies a bias correction to the abundance estimations obtained by closed 
  population models.
}
\usage{
closedp.bc(X, dfreq=FALSE, dtype=c("hist","nbcap"), t=NULL, t0=t,
           m=c("M0","Mt","Mh","Mth","Mb","Mbh"), h=NULL, h.control=list(), \dots)

\method{print}{closedp.bc}(x, \dots)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} for a description of the accepted formats).  }
  \item{dfreq}{ A logical. By default \code{FALSE}, which means that \code{X} has one row per unit. If TRUE, it indicates that the 
                matrix \code{X} contains frequencies in its last column.}
  \item{dtype}{ A characters string, either "hist" or "nbcap", to specify the type of data. "hist", the default, 
                means that \code{X} contains complete observed capture histories. "nbcap" means that \code{X} contains numbers 
                of captures (see \code{\link{Rcapture-package}} for details on data formats). If \code{m} is 
                "Mt", "Mth", "Mb" or "Mbh", \code{dtype} must be "hist". }
  \item{t}{ Requested only if \code{dtype="nbcap"}. A numeric specifying the total number of capture occasions in the experiment.}
  \item{t0}{ A numeric used for model M0 or for an Mh model other than Chao's lower bound model : 
             Models are fitted considering only the frequencies of units captured 
             1 to \code{t0} times. By default \code{t0=t}. }
  \item{m}{ A character string indicating the model to fit, either "M0"=M0 model, "Mt"=Mt model, "Mh"=Mh model, "Mth"=Mth model, 
            "Mb"=Mb model, "Mbh"=Mbh model.}
  \item{h}{ A character string ("LB", "Chao", "Poisson", "Darroch" or "Gamma") or a 
            numerical \code{R} function specifying the form of the column(s) for heterogeneity 
            in the design matrix.\cr
            "LB" and "Chao" represents Chao's lower bound model (the default),\cr 
            "Poisson" represents the function \eqn{f(k)=theta^k-1},\cr 
            "Darroch" represents the function \eqn{f(k)=k^2/2},\cr 
            and "Gamma" represents the function \eqn{f(k)=-\log(theta + k) + \log(theta)}{f(k)=-log(theta + k) + log(theta)}\cr
            where \eqn{k} is the number of captures and \eqn{theta} is a parameter specified with 
            the argument \code{theta}.\cr 
            If an \code{R} function is given, it must be the implementation of any convex 
            mathematical function \eqn{f(k)} (it must have only one argument).}
  \item{h.control }{ A list of elements to control the heterogeneous part of the model, if any.\cr
  For a Poisson or Gamma heterogeneous model:
\describe{
  \item{\code{theta}}{  The value of the parameter \eqn{theta} in \eqn{f(k)=theta^k-1} for the Poisson model
                (the default value is 2) and in \eqn{f(k)=-\log(theta + k) + \log(theta)}{f(k)=-log(theta + k) + log(theta)} 
                for the Gamma model (the default value is 3.5). }
}
For the Chao's lower bound Mth model:
\describe{
  \item{\code{neg}}{  If this option is set to \code{TRUE} (the default), negative eta parameters in 
              Chao's lower bound models are set to zero. }
}
}
  \item{\dots}{ Further arguments to be passed to \code{\link{glm}} or \code{\link{print.default}}.}
  \item{x}{ An object, produced by the \code{closedp.bc} function, to print.}
  }
\details{
For the Mt model:\cr 
When t=2, \code{closedp.bc} returns the Petersen estimator with Chapman's (1951) bias correction 
and the bias corrected standard error estimator of Seber (1970) and Wittes (1972).  
For t>2, \code{closedp.bc} implements the bias correction of Rivest and Levesque (2001).  
The estimate for N and its variance are calculated by solving an estimating equation as proposed 
in Seber (1982), not by fitting a Poisson regression.  This approach works for large values of t.
	
For other models:\cr
The bias correction is done through frequency modifications in Poisson regression as described 
in Rivest and Levesque (2001). The variances calculated with the modified frequencies are 
less biased than the standard ones, but they can overestimate the mean squared errors, 
especially when the data is sparse.

This function works with fairly large data set, except if an "Mth" model is requested. 
In this case, only heterogeneity of the form "LB", "Chao", "Poisson" with \code{theta=2} or "Darroch" 
is accepted.
}
\value{
\item{n}{ The number of captured units}
\item{t}{ The total number of capture occasions in the data matrix \code{X}.}
\item{t0}{ For models M0 and Mh only: the value of the argument \code{t0} used in the computations.}
\item{results}{ A table containing, for the fitted model:
    \describe{
    \item{\code{abundance}:}{  the estimated population size,}
    \item{\code{stderr}:}{  the standard error of the estimated population size,}
    \item{\code{infoFit}:}{  a numerical code giving information about error or warnings 
          encountered when fitting the model (see \code{\link{Rcapture-package}} for details).}
    }}
\item{glm.warn}{ Only if the corrected population size estimation was obtained with \code{\link{glm}}:
                 A vector of character strings. If the \code{\link{glm}} function generates
                 one or more warnings when fitting the model, a copy of these warnings are
                 stored in \code{glm.warn}. \code{NULL} if \code{\link{glm}} did not produce
                 any warnings. }
\item{neg.eta}{ For Chao's lower bound model Mth only: the position of the eta parameters set 
                to zero in the loglinear parameter vector, if any. } 
}
\note{ This function uses the \code{\link{glm}} function of the \pkg{stats} package, except for models Mt and Mh Chao's lower bound for which exact calculation is performed. }
\references{ 
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \doi{10.18637/jss.v019.i05}.

Chapman, D. G. (1951) Some properties of the hypergeometric distribution with applications to zoological sample censuses. \emph{University of California Publications in Statistics}, \bold{1}(7), 131-160.

Rivest, L.P. and Levesque, T. (2001) Improved loglinear model estimators of abundance in capture-recapture experiments. \emph{Canadian Journal of Statistics}, \bold{29}, 555-572.

Seber, G.A.F. (1970) The effects of trap response on tag recapture estimates. \emph{Biometrics}, \bold{26}, 13-22.

Seber, G.A.F. (1982) \emph{The Estimation of Animal Abundance and Related Parameters, 2nd edition}. New York: Macmillan.

Wittes, J.T. (1972) On the bias and estimated variance of Chapman's two-sample capture-recapture population estimate. \emph{Biometrics}, \bold{28}, 592-597.
}
\author{ Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and Sophie Baillargeon }

\seealso{\code{\link{closedp}} }
\examples{
# Third primary period of mvole data set
period3 <- mvole[, 11:15]
closedp.bc(period3, m = "Mh", h = "Darroch")
closedp.bc(period3, m = "Mh", h = "Gamma", h.control = list(theta = 3.5))

# BBS2001 data set
closedp.bc(BBS2001, dfreq = TRUE, dtype = "nbcap", t = 50, t0 = 20,
           m = "Mh", h = "Gamma", h.control = list(theta = 3.5))

# Seber (1982) p.107
# When there are 2 capture occasions, only models M0 and Mt can be fitted
X <- matrix(c(1,1,167,1,0,781,0,1,254), byrow = TRUE, ncol = 3)
closedp.bc(X, dfreq = TRUE, m = "M0")
closedp.bc(X, dfreq = TRUE, m = "Mt") 
}
\keyword{ models }
