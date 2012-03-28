\name{Rcapture-package}
\alias{Rcapture-package}
\alias{Rcapture}
\docType{package}
\title{ Loglinear Models for Capture-Recapture Experiments }
\description{
   Estimation of abundance and other demographic parameters for closed populations, open populations and the robust design in capture-recapture experiments using loglinear models.
}
\details{
\tabular{ll}{
Package: \tab Rcapture\cr
Type: \tab Package\cr
Version: \tab 1.2-1\cr
Date: \tab 2012-03-28\cr
License: \tab GPL-2\cr
}

DESCRIPTION OF THE ACCEPTED DATA SET FORMATS

A capture-recapture data set is given to the various \pkg{Rcapture} functions through the \code{X} argument. \code{X} must be a numeric matrix. The arguments \code{dfreq} and \code{dtype} indicate the format of the matrix. Each have two possible values, meaning that a total of four data set formats are possible with \pkg{Rcapture}.

1- If \code{dfreq=FALSE} and \code{dtype="hist"}, the default, \code{X} has one row per unit captured in the experiment. Each row is an observed capture history. It must contain only zeros and ones; the number one indicates a capture. In this case, the number of columns in the table represents the number of capture occasions in the experiment (noted \eqn{t}). Here is a fictive example of a data set of this type for \eqn{t=2}:\cr
1 1\cr
1 1\cr
1 0\cr
1 0\cr
1 0\cr
1 0\cr
0 1\cr

2- If \code{dfreq=TRUE} and \code{dtype="hist"}, \code{X} contains one row per capture history followed by its frequency. In that case, \code{X} has \eqn{t}+1 columns. As for the format presented previously, the first \eqn{t} columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. In this format, the example data set is represented by the following matrix:\cr
1 1 2\cr
1 0 4\cr
0 1 1\cr  

3- If \code{dfreq=FALSE} and \code{dtype="nbcap"}, \code{X} is simply a vector of numbers of captures. The length of the vector is \eqn{n}, the number of captured units. In this format, the example data set looks like:\cr
2 2 1 1 1 1 1\cr

4- If \code{dfreq=TRUE} and \code{dtype="nbcap"}, \code{X} is a 2 columns matrix. The first column contains the numbers of captures, the second columns contains the observed frequencies. In this format, the example data is:\cr
2 2\cr
1 5\cr

Only few functions have the \code{dtype} argument. Functions without \code{dtype} argument accept only a data matrix \code{X} of the form \code{dtype="hist"}. So the first two formats listed above are the most commun.

Formats with \code{dtype="nbcap"} are useful for experiments with a large number of capture occasions \eqn{t}. Often, no units will be caugth a large number of times, and the data set will contain no observations for \eqn{t} captures. Therefore, the number of capture occasions \eqn{t} cannot be deduced from \code{X} as it can be when \code{dtype="hist"}. So if one gives a data matrix \code{X} with \code{dtype="nbcap"}, one must also provide \code{t}, the number of capture occasions, as an additional argument.  

For now, the data formats with \code{dtype="nbcap"} are not generalized to the robust design. So \code{dtype} is not an argument of the \code{robustd.0} function.
}
\author{ Sophie Baillargeon \email{Sophie.Baillargeon@mat.ulaval.ca} and\cr
Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} }
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.

Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics}, \bold{43}(4), 783--791.

Cormack, R. M. (1985) Example of the use of glim to analyze capture-recapture studies. In \emph{Lecture Notes in Statistics 29: Statistics in Ornithology}, Morgan, B. et North, P. editors, New York,: Springer-Verlag, 242--274.

Cormack, R. M. (1989) Log-linear models for capture-recapture. \emph{Biometrics}, \bold{45}, 395--413.

Cormack, R. M. (1992) Interval estimation for mark-recapture studies of closed populations. \emph{Biometrics}, \bold{48}, 567--576.

Cormack, R. M. (1993) Variances of mark-recapture estimates. \emph{Biometrics}, \bold{49}, 1188--1193.

Cormack, R. M. and Jupp, P. E. (1991) Inference for Poisson and multinomial models for capture-recapture experiments. \emph{Biometrika}, \bold{78}(4), 911--916.

Rivest, L.P. and Levesque, T. (2001) Improved log-linear model estimators of abundance in capture-recapture experiments. \emph{Canadian Journal of Statistics}, \bold{29}, 555--572.

Rivest, L.P. and Daigle, G. (2004) Loglinear models for the robust design in mark-recapture experiments. \emph{Biometrics}, \bold{60}, 100--107.

Rivest, L.P. and Baillargeon, S. (2007) Applications and extensions of Chao's moment estimator for the size of a closed population. \emph{Biometrics}, \bold{63}(4), 999--1006.

Rivest, L.P. (2008) Why a time effect often has a limited impact on capture-recapture estimates in closed populations. \emph{Canadian Journal of Statistics}, \bold{36}(1), 75--84. 
}
\keyword{ package }
