\name{bunting}
\alias{bunting}
\docType{data}
\title{  Lazuli Bunting Data }
\description{
  This data set contains frequencies of capture histories for 8 years of observation of lazuli buntings. 
}
\usage{bunting}
\format{
255 by 9 numeric matrix, with the following columns:
  \describe{
    \item{\code{p1}, \code{p2}, \code{p3}, \code{p4}, \code{p5}, \code{p6}, \code{p7}, \code{p8}}{ Capture histories from eight periods}
    \item{\code{freq}}{ Observed frequencies for each capture history}
    }
}
\details{
  The data come from an eight-year (1973 to 1980) study by Allen W. Stokes of lazuli buntings wintering in Logan, Utah. 
  This data set is analyzed in Burnham and al. (1987) and in Cormack (1993).

  Each row of this data set represents an observed capture history followed by its frequency. 
}
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \doi{10.18637/jss.v019.i05}.

Burnham, K.P., Anderson, D.R., White, G.C., Brownie, C. and Pollock, K.H. (1987) Design and analysis methods for fish survival experiments based on release-recapture. \emph{American Fisheries Society Monographs 5}. Bethesda, Maryland.

Cormack, R. M. (1993) The flexibility of GLIM analyses of multiple recapture or resighting data. In J.D. Lebreton and P. North editors: \emph{Marked Individuals in the Study of Bird Population}. Basel, Switzerland: Birkhauser Verlag, 39--49. 
}
\examples{
descriptive(bunting, dfreq = TRUE)
  # 1430 birds among the 1681 birds seen (85\%) were caught only once. 
  # This suggests the presence of transient birds at each capture occasion. 

op.m1 <- openp(bunting, dfreq = TRUE)
op.m1$model.fit
plot(op.m1)
  # The residuals plot shows large residuals for the birds caught twice or
  # more while the residuals are small for birds caught once. The Jolly-Seber 
  # model does not fit well and the likely presence of transients might 
  # cause that. Let's remove the birds caught only once from the analysis.

keep2 <- rowSums(histpos.t(8)) > 1
op.m2<- openp(bunting, dfreq = TRUE, keep = keep2)
op.m2$model.fit
  # The deviance drop of 94 for 6 degrees of freedom is highly significant.
plot(op.m2)
  # The residual plot still shows Pearson residuals larger than 4. We redo
  # the analysis without the transient birds and without the large residuals.

keep3p <- residuals(op.m2$glm, type = "pearson") < 4
num3 <- ((1:255)[keep2])[keep3p]
keep3 <- rep(FALSE, 255)
keep3[num3] <- TRUE
op.m3 <- openp(bunting, dfreq = TRUE, keep = keep3)
cbind(op.m2$survivals, op.m3$survivals)
  # These changes have little impact on the survival estimates. 

  # We now test the equality of the survival probabilities and estimate its 
  # common value phi. Least squares estimates of phi and its standard error:
siginv <- solve(op.m2$cov[8:12, 8:12])
phi <- t(rep(1, 5)) \%*\% siginv \%*\% op.m2$survivals[2:6, 1] / 
      (t(rep(1, 5)) \%*\% siginv \%*\% rep(1, 5))
se <- 1 / sqrt(t(rep(1, 5)) \%*\% siginv \%*\% rep(1, 5))
cbind(estimate = phi, stderr = se)
  # The chi-square goodness of fit statistic for a constant survival 
  # and its pvalue are:
chisq4 <- t(op.m2$survivals[2:6, 1] - as.vector(phi)) \%*\% siginv \%*\% 
          (op.m2$survivals[2:6, 1] - as.vector(phi))
cbind(chi2stat = chisq4, pvalue = 1 - pchisq(chisq4, df = 4))
  # The hypothesis of a constant survival is accepted.
}
\keyword{datasets}
