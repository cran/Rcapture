\name{Rcapture-package}
\alias{Rcapture-package}
\alias{Rcapture}
\docType{package}
\title{ Loglinear Models for Capture-Recapture Experiments }
\description{
   Estimation of abundance and of other demographic parameters for closed populations, open populations and the robust design in capture-recapture experiments using loglinear models.
   
   This package focuses on closed populations. Since version 1.2-0, no new features have been added to open populations and robust design functions.
}
\details{
\tabular{ll}{
Package: \tab \pkg{Rcapture}\cr
Type: \tab Package\cr
Version: \tab 1.4-2\cr
Date: \tab 2014-08-28\cr
License: \tab GPL-2\cr
}

SUMMARY OF \pkg{Rcapture} CONTENTS

The \pkg{Rcapture} package contains nine capture-recapture data sets and the following functions:
\itemize{
\item Model fitting functions for\cr
> closed populations:
\itemize{
\item \code{\link{closedp}} functions: fit various loglinear models for abundance estimation, 
\item \code{\link{closedpCI}} functions: fit one customized loglinear model and calculate a confidance interval for the abundance estimation,
\item \code{\link{closedpMS.t}}: fits various hierarchical loglinear models in a perspective of model selection, 
\item \code{\link{closedp.bc}}: performs bias corrections to the abundance estimations from customized loglinear models,
\item \code{\link{closedp.Mtb}}: fits model Mtb, which cannot be fitted by any other function, for abundance estimation;
}
> open populations:
\itemize{
\item \code{\link{openp}}: computes various demographic parameters using a loglinear model;
}
> robust design: 
\itemize{
\item \code{\link{robustd}} functions: compute various demographic parameters and capture probabilities per period using a loglinear model.
}
\item Descriptive statistics functions: 
\itemize{
\item \code{\link{descriptive}}: produces basic descriptive statistics for capture-recapture data;
\item \code{\link{uifit}}: produces fit statistics concerning the \eqn{u_i}{ui}, i.e. the numbers of first captures on each capture occasion, for closed population models.
}
\item Data manipulation functions: 
\itemize{
\item \code{\link{histpos}} functions: builds a matrix of observable capture histories;
\item \code{\link{periodhist}}: merges capture occasions.\cr
}
}

DESCRIPTION OF DATA SET FORMATS

In capture-recapture experiments, the data collected consist of capture histories for captured units. A capture history is simply a serie of capture indicators for each capture event in the experiment. The capture history of one unit is expressed as a length \eqn{t} vector \eqn{w = (w_1, \ldots, w_t)}{w = (w1, ..., wt)}, where \eqn{w_j = 1}{wj = 1} if the unit is captured at the jth occasion and 0 if not. For closed populations, capture events are named capture occasions, whereas they are named capture periods for open populations.

Capture-recapture data sets are given to \pkg{Rcapture} functions through the \code{X} argument. \code{X} must be a numeric matrix. Arguments \code{dfreq} and \code{dtype} indicate the format of the matrix. Each have two possible values, meaning that four data set formats are possible with \pkg{Rcapture}.\cr

FORMAT 1 - CAPTURE HISTORY PER UNIT\cr
If \code{dfreq=FALSE} and \code{dtype="hist"} (the default), \code{X} has one row per unit captured in the experiment. Each row is an observed capture history. It must contain only zeros and ones; the number one indicates a capture. In this case, the number of columns in the table represents the number of capture occasions in the experiment (noted \eqn{t}). Here is an example of a data set of this type for \eqn{t=2}:\cr
1 1\cr
1 1\cr
1 0\cr
1 0\cr
1 0\cr
1 0\cr
0 1\cr

FORMAT 2 - AGGREGATED CAPTURE HISTORIES\cr
If \code{dfreq=TRUE} and \code{dtype="hist"}, \code{X} contains one row per observed capture history followed by its frequency. In that case, \code{X} has \eqn{t}+1 columns. As for format 1, the first \eqn{t} columns of \code{X}, identifying the capture histories, must contain only zeros and ones. The number one indicates a capture. In this format, 
the example data set is represented by the following matrix:\cr
1 1 2\cr
1 0 4\cr
0 1 1\cr
If a possible capture history is not observed, it can appear in \code{X} with a frequency of zero, or it can simply be omitted.\cr

FORMAT 3 - NUMBER OF CAPTURES PER UNIT\cr
If \code{dfreq=FALSE} and \code{dtype="nbcap"}, \code{X} is a vector with the number of captures for every captured unit. Therefore, this format does not contain complete capture histories. Instead, capture histories are summarized through the number of captures. In this format, the example data set looks like:\cr
2 2 1 1 1 1 1\cr

FORMAT 4 - AGGREGATED NUMBERS OF CAPTURES\cr
If \code{dfreq=TRUE} and \code{dtype="nbcap"}, \code{X} is a 2 columns matrix. The first column contains the observed numbers of captures, the second columns contains their frequencies. 
In this format, the example data is:\cr
2 2\cr
1 5\cr


DETAILS ABOUT FORMATS WITH NUMBERS OF CAPTURES

Only few functions have the \code{dtype} argument. Functions without \code{dtype} argument accept only a data matrix \code{X} of the form \code{dtype="hist"}. So the first two formats listed above are the most common.

Formats with \code{dtype="nbcap"} are used for captures in continuous time (see below). They are also useful to reduce the size of the data set for experiments with a large number of capture occasions \eqn{t} (often with no units caught a large number of times). For theses formats, the number of capture occasions \eqn{t} cannot be deduced from \code{X} as it can be with \code{dtype="hist"}. One has no garanties that the larger number of captures observed is the total number of capture occasions. Therefore, if one gives a data matrix \code{X} with \code{dtype="nbcap"}, one must also provide \code{t}, the number of capture occasions, as an additional argument.  

For now, the data formats with \code{dtype="nbcap"} are not generalized to the robust design. So \code{dtype} is not an argument of the \code{robustd.0} function.\cr


CAPTURES IN CONTINUOUS TIME

In some capture-recapture experiments, there is no well defined capture occasions. 
Captures occur in continuous time. The data set \code{\link{ill}} comes from
such an experiment. Bohning and Schon (2005) call this type of capture-recapture 
data \emph{repeated counting data}. These data sets always have the format \code{dtype="nbcap"}. 

We can estimate abundance for data of this type using the option \code{t=Inf} with the
functions \code{\link{closedpCI.0}} and \code{\link{closedpCI.0}}. The function 
\code{\link{descriptive}} also accepts \code{t=Inf}. It modifies the y coordinate
of the exploratory heterogeneity graph.\cr


DISTINCTION BETWEEN \code{.t} and \code{.0} FUNCTIONS

Capture recapture models for closed population aim at estimating the population size by modelling the probabilities of the different capture histories. The data available to fit these models consist of observed frequencies of capture histories. These frequencies are modeled in Rcapture using loglinear models for frequency tables. For functions with a name ending with \code{.t} (\code{\link{closedp.t}}, \code{\link{closedpCI.t}}, \code{\link{closedpMS.t}} and \code{\link{robustd.t}}), the observed values of the response variable is the vector of frequencies for every observable capture history. It has length \eqn{2^t-1}.

However, for a model without temporal effect (assuming that the probability of capturing a unit does not vary between capture occasions), all the information needed to fit the model can be found in aggregated data. Functions with a name ending with \code{.0} (\code{\link{closedp.0}}, \code{\link{closedpCI.0}} and \code{\link{robustd.0}}) fit models using as response variable the number of units captured \eqn{i} times, for \eqn{i=1,\ldots,t}. It is a vector of length \eqn{t}, which is much shorter than \eqn{2^t-1} for a large \eqn{t}. Because of an appropriate offset added to the model, the results from a \code{.0} function match exaclty the results from a \code{.t} function for a corresponding capture-recapture model. Because \code{.0} functions deal with smaller design matrix, they run faster than \code{.t} functions, but they cannot fit models with a temporal effect. \cr


FUNCTIONS USED FOR MODEL FITTING

Most of the \pkg{Rcapture} functions use the function \code{\link{glm}} of the \pkg{stats} package to fit a loglinear model. However, the function \code{\link{optim}}, again of the \pkg{stats}, is also used by three functions:
\itemize{
\item \code{\link{closedpCI.t}}, \code{\link{closedpCI.0}} and \code{\link{closedpMS.t}}: when a normal heterogeneous models is requested,
\item \code{\link{closedp.Mtb}}: because model Mtb does not have a loglinear form.\cr 
}


ERRORS AND WARNINGS MANAGEMENT IN CLOSED POPULATION FUNCTIONS

If an error occurs while executing a function fitting only one closed population model (\code{\link{closedpCI}} functions, \code{\link{closedp.bc}}, \code{\link{closedp.Mtb}}), the execution is stopped and the error message is printed (usual behavior in R). However, if an error occurs while fitting a model in a call to a closedp population function fitting more than one model (\code{\link{closedp}} functions, \code{\link{closedpMS.t}}), the execution of the call is not stopped. Instead, the row in the results table for the problematic model is filled with \code{NA} and the error message is stored in an output value. This value is called \code{glm.err} for \code{\link{closedp}} functions. It is however called \code{fit.err} for \code{\link{closedpMS.t}} because this function do not always use \code{\link{glm}} to fit the model (as mentionned above, \code{\link{optim}} is used for normal heterogeneous models).

Warning messages while fitting a closed population model, if any, are stored in an output value called 'glm.warn', 'fit.warn' or 'optim.warn' depending on the function. They are not printed in the console. To inform the user that a warning occured, the last column of the results table, named \code{infoFit}, contains a numerical code giving information about errors or warnings encountered. Here is a description of the meaning of the numerical code:
\describe{
\item{   0}{  no error or warning occured while fitting the model;}
\item{  -1}{  an error occured while fitting the model;}
\item{   1}{  a warning indicating that the model fit is questionnable occured (algorithm did not converge, non-positive sigma estimate for a normal heterogeneous model or large asymptotic bias);} 
\item{   2}{  the warning 'design matrix not of full rank' occured, therefore some model's coefficients are not estimable;}
\item{   3}{  a warning not of type 1 or 2 occured (the \code{\link{glm}} warning 'fitted rates nummerically 0 occured' is often encountered with small frequencies, it does not always mean that the model fit is questionnable).}
}
The elements in the column \code{infoFit} can contain more than one number since more than one warning can occur. For exemple, if \code{infoFit} takes the value \code{13} for a model, it means that at least one warning of type 1 and one warning of type 3 have occured. If more than one warning of the same type are encountered, the number representing the type of warning is not repeated in \code{infoFit}. 

}
\author{ Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and Sophie Baillargeon }

\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.

Bohning, D. and Schon, D. (2005) Nonparametric Maximum Likelihood Estimation of Population Size Based on the Counting Distribution. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \bold{54}(4), 721-737.

Chao, A. (1987) Estimating the population size for capture-recapture data with unequal catchability. \emph{Biometrics}, \bold{43}(4), 783--791.

Cormack, R. M. (1985) Example of the use of glim to analyze capture-recapture studies. In \emph{Lecture Notes in Statistics 29: Statistics in Ornithology}, Morgan, B. et North, P. editors, New York,: Springer-Verlag, 242--274.

Cormack, R. M. (1989) Loglinear models for capture-recapture. \emph{Biometrics}, \bold{45}, 395--413.

Cormack, R. M. (1992) Interval estimation for mark-recapture studies of closed populations. \emph{Biometrics}, \bold{48}, 567--576.

Cormack, R. M. (1993) Variances of mark-recapture estimates. \emph{Biometrics}, \bold{49}, 1188--1193.

Cormack, R. M. and Jupp, P. E. (1991) Inference for Poisson and multinomial models for capture-recapture experiments. \emph{Biometrika}, \bold{78}(4), 911--916.

Rivest, L.P. and Levesque, T. (2001) Improved loglinear model estimators of abundance in capture-recapture experiments. \emph{Canadian Journal of Statistics}, \bold{29}, 555--572.

Rivest, L.P. and Daigle, G. (2004) Loglinear models for the robust design in mark-recapture experiments. \emph{Biometrics}, \bold{60}, 100--107.

Rivest, L.P. and Baillargeon, S. (2007) Applications and extensions of Chao's moment estimator for the size of a closed population. \emph{Biometrics}, \bold{63}(4), 999--1006.

Rivest, L.P. (2008) Why a time effect often has a limited impact on capture-recapture estimates in closed populations. \emph{Canadian Journal of Statistics}, \bold{36}(1), 75--84. 
}

\examples{
# Here is an example on the lesbian data set.

data(lesbian)

desc <- descriptive(lesbian, dfreq = TRUE)
desc
plot(desc)

# 1612 out of 2185 individuals (74\%) appear on one list only.
# The exploratory heterogeneity graph are not quite linear.
# Some heterogeneity in the units capture probabilities
# seem present in the data set.

closedp(lesbian, dfreq = TRUE)

# According to the BIC, the best model is Mth Darroch.
# Let's see if adding  interactions between capture 
# histories to the model could improve the model's fit.

closedpMS.t(lesbian, dfreq = TRUE, h = "Darroch")

# According to the BIC, the best heterogeneous Darroch model
# contains the double interactions 12, 13, 14. 
# Here is the profile likelihood confidence interval for the
# abundance estimation from this model.

closedpCI.t(lesbian, dfreq = TRUE, mX = "[12,13,14]", h = "Darroch")


####################################################

# Example to illustrate warnings management in closed population functions.

# Here is a capture-recapture data set one could encounter.

crdata <- cbind(histpos.t(4), c(0,0,3,0,0,0,0,0,0,0,0,1,0,0,2))

# This data set contains 4 capture occasions but only 6 captured units.
# Fitting capture-recapture models on this data set is quite useless.
# The population size should be very close to the sample size.

# Such small frequencies in a capture-recapture data set should
# lead to warnings when fitting a loglinear model on it.

ex <- closedp.t(crdata, dfreq=TRUE)
ex

# Many models produce warnings of type 1 indicating that the model fit 
# is questionnable. The very large abundance estimation for some models 
# are another indicator of questionable model fits.
# Details about the warnings are found in the glm.warn element of the output.

ex$glm.warn
}

\keyword{ package }
