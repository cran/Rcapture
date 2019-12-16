\name{ill}
\alias{ill}
\docType{data}
\title{
Illegal immigrants in the Netherlands
}
\description{
This data set contains the frequency \eqn{f_i}{fi} that an illegal immigrant who is not
effectively expelled is apprehended exactly \eqn{i} times.
}
\usage{ill}
\format{
6 by 2 numeric matrix, with the following columns:
  \describe{
    \item{\code{nbcap}}{ Numbers of captures, i.e. number of times an illegal immigrant is apprehended by the police}
    \item{\code{freq}}{ Observed frequencies for each number of captures}
    }
}
\details{
Bohning and Schon (2005) presented this data set as follows:\cr
The illegal immigrant data are from van der Heijden et al. (2003) and refer to count data on
illegal immigrants in four large cities in the Netherlands who could not be effectively expelled
from the country. Illegal immigrants who are apprehended by the police often cannot be effectively
expelled because either they refuse to disclose their nationality or their home country does
not cooperate in receiving them back. In such cases they will be asked by the police to leave
the country, but it is unlikely that they will abide by this request. The data were collected by the
police and date back to the year 1995.
}
\source{
van der Heijden, P.G.M., Cruyff, M. and van Houwelingen, H.C. (2003) Estimating the size of a criminal population from police records using the truncated Poisson regression model. \emph{Statistica Neerlandica}, \bold{57}(3), 289-304.
}
\references{
Bohning, D. and Schon, D. (2005) Nonparametric Maximum Likelihood Estimation of Population Size Based on the Counting Distribution. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \bold{54}(4), 721-737.
}
\examples{
sdesc <- descriptive(ill, dtype = "nbcap", dfreq = TRUE, t = Inf)
plot(sdesc)
   # A mixture model looks appropriate

closedp.0(ill, dtype = "nbcap", dfreq = TRUE, t = Inf)

   # We can try to fit a normal mixture model:
closedpCI.0(ill, dtype = "nbcap", dfreq = TRUE, t = Inf, m = "Mh", h = "Normal")
   # We get an estimate similar to the Mh Gamma3.5 estimate.

   # Estimates are highly variable and it seems difficult to come up with a
   # definitive answer. The lower bound estimate is useful in this context.
closedpCI.0(ill, dtype = "nbcap", dfreq = TRUE, t = Inf, m = "Mh", h = "LB")
   # Considering the lower limit of a 95\% confidence estimate for the lower bound, 
   # there should be at least 8 000 illegal immigrants in the Netherlands. 
   # Less than 25\% have been caught.
}
\keyword{datasets}
