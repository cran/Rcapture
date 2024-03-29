\name{uifit}
\alias{uifit}
\title{ Closed Population Models Fit Statistics Concerning the ui }
\description{
  This function produces fit statistics concerning the \eqn{u_i}{ui}, i.e. the numbers of first captures on each capture occasion, 
  for closed population models. It also forecasts, for some models, the number of additional units that would be captured if the 
  experiment was continued for five more occasions.
}
\usage{
uifit(x.closedp)
}
\arguments{
  \item{x.closedp}{ An object produced by the \code{closedp.t} function. }
}
\value{
  \item{predicted }{ The observed and predicted values of the ui frequencies, \eqn{i = 1,\ldots,t}, for each model in \code{closedp}. 
                     For models M0, Mh Poisson2, Mh Darroch and Mb, predicted values are also extrapolated for \eqn{i = t+1,\ldots,t+5}.}
  \item{fit.stat }{ Chi-square fit statistics for each model in \code{closedp} 
  \eqn{(\sum_{i=1}^t \frac{(\mbox{observed}_i -\mbox{predicted}_i)^2}{\mbox{predicted}_i})}{(sum((observed - predicted)^2/predicted))}.}
  \item{day.first.capt }{ The mean and variance of the day of first capture, calculated with the observed and predicted \eqn{u_i}{ui}, 
                          for each model in \code{closedp.t}. }
}
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \doi{10.18637/jss.v019.i05}.
}
\author{ Louis-Paul Rivest \email{Louis-Paul.Rivest@mat.ulaval.ca} and Sophie Baillargeon }

\seealso{\code{\link{closedp}} }
\examples{
# Third primary period of mvole data set
period3 <- mvole[, 11:15]
cp <- closedp(period3)
uifit(cp)
}
\keyword{ univar }
