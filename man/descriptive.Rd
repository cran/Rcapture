\name{descriptive}
\alias{descriptive}
\alias{print.descriptive}
\alias{plot.descriptive}
\title{ Descriptive Statistics }
\description{
  This function produces descriptive statistics for capture-recapture data.
}
\usage{
descriptive(X, dfreq=FALSE)

\method{print}{descriptive}(x, \dots)

\method{plot}{descriptive}(x, \dots)
}
\arguments{
  \item{X}{ The table of the observed capture histories in one of the two accepted formats. In the default format, it has one row per unit captured in the experiment. In this case, the number of columns in the table represents the number of capture occasions in the experiment (noted \eqn{t}). In the alternative format, it contains one row per capture history followed by its frequency. In that case, \code{X} has \eqn{t}+1 columns. The first \eqn{t} columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. }
  \item{dfreq}{ This argument specifies the format of the data matrix \code{X}. By default, it is set to FALSE, which means that \code{X} has one row per unit. If it is set to TRUE, then the matrix \code{X} contains frequencies in its last column.}
  \item{x}{ An object, produced by the \code{descriptive} function, to print or to plot.}
  \item{\dots}{ Further arguments passed to or from other methods.}
}
\details{ 
  The \code{plot.descriptive} function produces exploratory heterogeneity graphs. In absence of heterogeneity, the relations presented in the two graphs should be almost linear. Convex fonctions indicate heterogeneity.  
}
\value{
  \item{n}{ The total number of captured units.}
  \item{base.freq }{ A table of basic descriptive statistics. For \eqn{i=1,\ldots,t} :
    \item{fi}{the number of units captured \eqn{i} times}
    \item{ui}{the number of units captured for the first time on occasion \eqn{i} }
    \item{vi}{the number of units captured for the last time on occasion \eqn{i} }
    \item{ni}{the number of units captured on occasion \eqn{i} }.
    }
  \item{m.array }{ A \eqn{t} by \eqn{t} matrix containing recapture statistics for units released on occasion \eqn{i}, for \eqn{i=1,\ldots,t} :
    \item{ni}{the number of units captured on occasion \eqn{i} }
    \item{c2}{the number of units caught on occasion \eqn{i} recaptured for the first time on occasion 2}
    \cr ....
    \item{ct}{the number of units caught on occasion \eqn{i} recaptured for the first time on occasion \eqn{t} }
    \item{not recapt}{the number of units caught on occasion \eqn{i} never recaptured}.
    }
}
\references{ 
Baillargeon, S. and Rivest, L.P. (2007). The Rcapture package: Loglinear models for capture-recapture in R. \emph{Journal of Statistical Software}, to appear (available online at \url{http://www.jstatsoft.org/}).

Lindsay, B. G. (1986) Exponential family mixture models (with least-squares estimators).  \emph{ Annals of Statistics}, 14, 124-137.

Rivest, L.P. (2007). Why a time effect has a limited impact on capture-recapture estimates in closed populations. \emph{The Canadian Journal of Statistics}, under revision.
}
\author{Sophie Baillargeon \email{sbaillar@mat.ulaval.ca} and Louis-Paul Rivest \email{lpr@mat.ulaval.ca} }
\examples{
data(hare)
desc<-descriptive(hare)
desc
plot(desc)
}
\keyword{ univar }
