\name{duck}
\alias{duck}
\docType{data}
\title{  Eider Duck Data }
\description{
  This data set contains capture-recapture data for eider ducks. 
}
\usage{duck}
\format{
63 by 7 numeric matrix, with the following columns:
  \describe{
    \item{\code{p1}, \code{p2}, \code{p3}, \code{p4}, \code{p5}, \code{p6}}{ Capture histories from six periods}
    \item{\code{freq}}{ Observed frequencies for each capture history}
    }
}
\details{
  The data points are extracted from a 25-year study by Coulson (1984). The capture periods are six consecutive years : years 19-24. 
  This data set is analyzed in Cormack (1989).
  
  Each row of this data set represents an observed capture history followed by its frequency. 
}
\source{
Coulson, J. C. (1984) The population dynamics of the Eider Duck \emph{Somateria mollissima} and evidence of extensive non breeding by adults ducks. \emph{Ibis}, \bold{126}, 525--543.
}
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \doi{10.18637/jss.v019.i05}.

Cormack, R. M. (1989) Loglinear models for capture-recapture. \emph{Biometrics}, \bold{45}, 395--413.
}
\examples{
op.m1 <- openp(duck, dfreq = TRUE)
op.m1$model.fit
  # The pvalue of the goodness of fit test based on the deviance is
1 - pchisq(op.m1$model.fit[1, 1], df = 49)
plot(op.m1)
  # The residual plot shows a large residual for the 13 individuals 
  # captured all the times. We redo the analysis without them.

keep2 <- rowSums(histpos.t(6)) != 6
op.m2 <- openp(duck, dfreq = TRUE, keep = keep2)
op.m2$model.fit
1 - pchisq(op.m2$model.fit[1, 1],df = 48)
  # The fit is still not satisfactory.
plot(op.m2)
  # The residual plot has the convex shape characteristic of 
  # heterogeneity in the capture probabilities. We also remove the 
  # individuals caught at 5 periods out of 6.

keep3 <- rowSums(histpos.t(6)) < 5
op.m3 <- openp(duck, dfreq = TRUE, keep = keep3)
op.m3$model.fit
1 - pchisq(op.m3$model.fit[1, 1], df = 42)
  # The fit is better but there is still heterogeneity in the data. 

  # To investigate whether the capture probabilities are homogeneous, 
  # one can fit a model with equal capture probabilities.
op.m4 <- openp(duck, dfreq = TRUE, m = "ep", keep = keep3)
op.m4$model.fit
  # It gives a much larger deviance; this model is not considered anymore.

  # We now investigate models for the growth rate N[i+1]/N[i] of this 
  # population using the multivariate normal distribution for the 
  # abundance estimates. The growth rates and their standard errors are
growth <- op.m3$N[3:5]/op.m3$N[2:4]
partial <- matrix(c(-op.m3$N[3]/op.m3$N[2]^2, 1/op.m3$N[2], 0, 0,
                    0, -op.m3$N[4]/op.m3$N[3]^2, 1/op.m3$N[3], 0,
		                0, 0, -op.m3$N[5]/op.m3$N[4]^2, 1/op.m3$N[4]), 3, 4, byrow = TRUE)
sig <- partial \%*\% op.m3$cov[9:12, 9:12] \%*\% t(partial)
cbind(estimate = growth, stderr = sqrt(diag(sig)))
  # An estimate for the common growth rate is
siginv <- solve(sig)
growth.e <- t(rep(1, 3)) \%*\% siginv \%*\% growth/(t(rep(1, 3)) \%*\% siginv \%*\% rep(1, 3))
se <- 1/sqrt(t(rep(1, 3)) \%*\% siginv \%*\% rep(1, 3))
cbind(estimate = growth.e, stderr = se)
  # A chi-square statistic for testing the equality of the growth rates 
  # and its pvalue
chisq2 <- t(growth - as.vector(growth.e)) \%*\% siginv \%*\% (growth - as.vector(growth.e))
c(chi2stat = chisq2, pvalue = 1 - pchisq(chisq2, df = 2))
  # The hypothesis of a common growth rate is rejected
}
\keyword{datasets}
