\name{cbird}
\alias{cbird}
\docType{data}
\title{
Catbird Site Occupancy Data
}
\description{
This data set contains site occupancy data for catbirds. The number of capture occasions \eqn{t} is 11.
}
\usage{cbird}
\format{
6 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{nbcap}}{ Numbers of captures}
    \item{\code{freq}}{ Observed frequencies for each number of captures}
    }
}
\details{
This data set is presented in Royle (2006). It comes from an experiment in which 50 locations have been visited 11 times. 
The presence or absence of catbirds was noted at each visit. The experiment aims at estimating the number of locations (<=50) 
at which catbirds are present. 

In the data set \code{cbird}, the frequencies represent the number of locations at which cbirdirbs have been detected \eqn{i} times, 
out of \eqn{t=11} visits. The frequencies for \eqn{i=7,\ldots,11} are null and do not appear in the data set. 
}
\source{
Royle, A. J. (2006) Site occupancy models with heterogeneous detection probabilities. \emph{Biometrics}, \bold{62}, 97--102 }
\examples{
desc <- descriptive(cbird, dfreq = TRUE, dtype = "nbcap", t = 11)
plot(desc)
   # The heterogeneity graph is mildly convex, mostly because of 
   # 2 sites where catbirds have been seen six times.

closedp.0(cbird, dfreq = TRUE, dtype = "nbcap", t = 11)
   # The residuals are OK. The occupancy rate (N/50) estimate 
   # varies between 40.6\% for M0 and 49.4\% for Mh Gamma.
   # One could also try fitting M0 removing the 2 unusual sites:
closedp.0(cbird, dfreq = TRUE, dtype = "nbcap", t = 11, t0 = 5)
   # M0 fits very well and the occupancy rate estimate is 42.6\%.
}
\keyword{datasets}
