\name{BBS2001}
\alias{BBS2001}
\docType{data}
\title{
Species Richness Data from the North American Breeding Bird Survey (BBS) in 2001
}
\description{
This data set contains species richness data from the North American Breeding Bird Survey (BBS) in 2001. 
The number of capture occasions \eqn{t} is 50.
}
\usage{BBS2001}
\format{
36 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{nbcap}}{ Numbers of captures}
    \item{\code{freq}}{ Observed frequencies for each number of captures}
    }
}
\details{
This data set is presented in Dorazio and Royle (2003). It comes from the North American Breeding Bird Survey (BBS). 
This survey includes more than 4000 39.4 km survey routes throughout North America. Along each survey route an observer 
stops at 50 equidistant locations and records the species identity and number of all birds heard or seen within a 3-minute period. 

The data in \code{BBS2001} was collected along Route 123, in Maryland, in 2001. In this data set a capture occasion is in 
fact a location and the frequencies represent the number of species detected at \eqn{i} locations out of the 50 locations 
on that route. The frequencies for \eqn{i=37,\ldots,50} are null and do not appear in the data set. 
}
\source{
Dorazio, R. M. and Royle, A. J. (2003) Mixture models for estimating the size of a closed population when capture rates vary among individuals. \emph{Biometrics}, \bold{59}, 351--364
}
\examples{
desc <- descriptive(BBS2001, dfreq = TRUE, dtype = "nbcap", t = 50)
plot(desc)
   # Note in this plot the convex shape typical of model Mh.

cp <- closedp.0(BBS2001, dfreq = TRUE, dtype = "nbcap", t = 50, t0 = 20)
cp
plot(cp)
   # The Gamma estimator has the smallest deviance and all 
   # its residuals are small.
   
   # Let's calculate a profile confidence interval for the gamma 
   # estimator:
closedpCI.0(BBS2001, dfreq = TRUE, dtype = "nbcap", t = 50, t0 = 20,
            m = "Mh", h = "Gamma", h.control = list(theta = 3.5))		

   # One can also calculate a profile confidence interval for Chao's 
   # lower bound:
closedpCI.0(BBS2001, dfreq = TRUE, dtype = "nbcap", t = 50, t0 = 20,
            m = "Mh", h = "Chao")   
}
\keyword{datasets}
