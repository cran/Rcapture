\name{descriptive}
\alias{descriptive}
\alias{print.descriptive}
\alias{plot.descriptive}
\title{ Descriptive Statistics }
\description{
  This function produces descriptive statistics for capture-recapture data.
}
\usage{
descriptive(X, dfreq=FALSE, dtype=c("hist","nbcap"), t=NULL)

\method{print}{descriptive}(x, \dots)

\method{plot}{descriptive}(x, main="Exploratory Heterogeneity Graph", \dots)
}
\arguments{
  \item{X}{ The matrix of the observed capture histories (see \code{\link{Rcapture-package}} for a description of the accepted formats).  }
  \item{dfreq}{ A logical. By default FALSE, which means that \code{X} has one row per unit. If TRUE, it indicates that 
                the matrix \code{X} contains frequencies in its last column.}
  \item{dtype}{ A characters string, either "hist" or "nbcap", to specify the type of data. "hist", the default, means that 
                \code{X} contains complete observed capture histories. "nbcap" means that \code{X} contains numbers of captures 
                (see \code{\link{Rcapture-package}} for details on data formats). }
  \item{t}{ Requested only if \code{dtype="nbcap"}. A numeric specifying the total number 
            of capture occasions in the experiment. The value \code{t=Inf}
            is accepted. It indicates that captures occur in continuous time.
            (see \code{\link{Rcapture-package}}).}
  \item{x}{ An object, produced by the \code{descriptive} function, to print or to plot.}
  \item{main}{ A main title for the plot}
  \item{\dots}{ Further arguments to be passed to methods (see \code{print.default} and \code{plot.default}).}
}
\details{ 
  The \code{plot.descriptive} function produces an exploratory heterogeneity graph. 
  In the absence of heterogeneity, the relation(s) presented in the graph should be almost linear. 
  Convex functions indicate heterogeneity (Baillargeon and Rivest (2007)).      
}
\value{
  \item{n}{ The total number of captured units.}
  \item{base.freq }{ A table of basic descriptive statistics. For \eqn{i=1,\ldots,t} :\cr
    fi: the number of units captured \eqn{i} times\cr
    ui: the number of units captured for the first time on occasion \eqn{i}\cr
    vi: the number of units captured for the last time on occasion \eqn{i}\cr
    ni: the number of units captured on occasion \eqn{i}\cr
    If \code{dtype="nbcap"}, only the the number of units captured \eqn{i} times, fi, can be computed.
    }
  \item{m.array }{ Only if \code{dtype="hist"}: A \eqn{t} by \eqn{t} matrix containing recapture statistics for units 
                   released on occasion \eqn{i}, for \eqn{i=1,\ldots,t} :\cr
    ni: the number of units released on occasion \eqn{i};\cr
    c2: the number of units released on occasion \eqn{i} and recaptured for the first time on occasion 2;\cr
    \dots\cr
    ct: the number of units released on occasion \eqn{i} and recaptured for the first time on occasion \eqn{t};\cr
    not recapt: the number of units released on occasion \eqn{i} and never recaptured again.\cr
    The m.array matrix is used in Mark (see White and Burnham 1999) and MSurge (see Choquet, Reboulet, Pradel, Gimenez, 
    and Lebreton 2004) to fit open population models to the data.
    }
   \item{call }{ The function call (object of class "call").}
}
\references{ 
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.

Choquet, R., Reboulet, A.M., Pradel, R., Gimenez, O. and Lebreton, J.D. (2004). M-Surge: New Software Specifically Designed for Multistate Capture-Recapture Models. \emph{Animal Biodiversity and Conservation}, \bold{27.1}, 207--215.

Hoaglin, D. C. (1980) A Poissonness Plot, \emph{The American Statistician}, \bold{34}, 146--149

Lindsay, B. G. (1986) Exponential family mixture models (with least-squares estimators).  \emph{ Annals of Statistics}, \bold{14}, 124--137.

Rivest, L.P. (2008) Why a time effect often has a limited impact on capture-recapture estimates in closed populations. \emph{Canadian Journal of Statistics}, \bold{36}(1), 75--84.

White, G. and Burnham, K.P. (1999) Program Mark: Survival Estimation from Populations of Marked Animals. \emph{Bird Study}, \bold{46} (Supplement), 120--138. The software can be downloaded here: \url{http://www.phidot.org/software/mark/}.
}
\author{ Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and Sophie Baillargeon }

\examples{
# hare data set
hare.desc <- descriptive(hare)
hare.desc
plot(hare.desc)

# cbird data set
cbird.desc1 <- descriptive(cbird, dfreq = TRUE, dtype = "nbcap", t = 11)
plot(cbird.desc1)
# To illustrate the option t = Inf.
cbird.desc2 <- descriptive(cbird, dfreq = TRUE, dtype = "nbcap", t = Inf)
plot(cbird.desc2)
# The y coordinate has changed.
}
\keyword{ univar }
