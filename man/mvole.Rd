\name{mvole}
\alias{mvole}
\docType{data}
\title{  Robust Design Data for Adult Male Meadow Voles }
\description{
  This data set contains robust design capture history data for adult male meadow voles. 
}
\usage{mvole}
\format{
171 by 30 numeric matrix, with the following columns:
  \describe{
    \item{\code{c11}, \code{c12}, \code{c13}, \code{c14}, \code{c15}}{Capture histories from the 5 capture occasions within primary period 1}
    \item{\code{c21}, \code{c22}, \code{c23}, \code{c24}, \code{c25}}{Capture histories from the 5 capture occasions within primary period 2}
    \item{\code{c31}, \code{c32}, \code{c33}, \code{c34}, \code{c35}}{Capture histories from the 5 capture occasions within primary period 3}
    \item{\code{c41}, \code{c42}, \code{c43}, \code{c44}, \code{c45}}{Capture histories from the 5 capture occasions within primary period 4}
    \item{\code{c51}, \code{c52}, \code{c53}, \code{c54}, \code{c55}}{Capture histories from the 5 capture occasions within primary period 5}
    \item{\code{c61}, \code{c62}, \code{c63}, \code{c64}, \code{c65}}{Capture histories from the 5 capture occasions within primary period 6}
  }
}
\details{
  The data set is extracted from Table 19.1 of Williams, Nichols and Conroy (2002). The capture occasions represent five 
  consecutive days of trapping each month from June to December 1981 at Patuxent Wildlife Research Center, Laurel, Maryland.

  Each row of this data set represents the capture history of one animal. 
}
\note{
  In this data set, ten animals are in fact not released after capture. These trap deaths are not identified. 
}
\source{
Williams, B.K., Nichols, J.D., and Conroy, M.J. (2002) \emph{Analysis and Management of Animal Populations}, San Diego: Academic Press.
}
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \url{http://www.jstatsoft.org/v19/i05}.
}
\examples{
  # This example deals only with the first four primary periods of the data set.
mvole4 <- mvole[, 1:20]

  # First, a between primary period Jolly-Seber analysis is obtained.
mvole4.pp <- periodhist(mvole4, vt = rep(5,4))
op.m1 <- openp(mvole4.pp, dfreq = TRUE)
plot(op.m1)
  # There is one large residual, removing the corresponding capture history 
  # from the analysis does not change the results.  The model fits well.
keep2 <- residuals(op.m1$glm, type = "pearson") < 4
op.m2 <- openp(mvole4.pp, dfreq = TRUE, keep = keep2)
op.m2$model.fit

  # To find a suitable model within each primary period, the function closedp.t
  # can be used repeatedly. Heterogeneity is detected in all periods except
  # the second one where the data collection was perturbed (the last capture
  # occasion doesn't have any new capture and is taken out of the analysis). 

  # In a robust design, we use Mh models for all primary periods bearing in
  # mind the questionable fit in the second one. Since there is no time effect 
  # within primary periods, we use the function robustd.0 to fit the model.
  ### The following command might take a few minutes to run.
rd.m1 <- robustd.0(mvole4[, -10], vt = c(5, 4, 5, 5), vm = "Mh", vh = "Chao")
rd.m1$model.fit
rd.m1$emig.fit
  # The test for temporary immigration is not significant meaning that capture 
  # probabilities estimated with the Jolly-Seber model are not different from 
  # those estimated with the individual closed population models. The 
  # differences, on the logit scale, of the Jolly-Seber minus the closed 
  # population models capture probabilities are
rd.m1$emig.param
  # Even in period 2, where the closed population model does not fit well, the 
  # difference on the logit scale is non significant (estimate=.56, s.e.=1.13). 

  # The following command allows to fit a robust design that does not specify 
  # any model for the second period.
  ### The following command might take a few minutes to run.
rd.m3 <- robustd.0(mvole4[, -10], vt = c(5, 4, 5, 5), 
                   vm = c("Mh", "none", "Mh", "Mh"), vh = "Chao")

  # With Darroch's model, the closed population estimates of the capture 
  # probabilities are significantly smaller than those obtained from the 
  # Jolly-Seber model.  This cannot be interpreted as indicating temporary 
  # immigration. This suggests that Darroch's model is not appropriate within
  # primary sessions.

  # The smallest AIC is obtained with the Poisson model, with parameter a=1.5 
  # within sessions.
rd.m4 <- robustd.0(mvole4[, -10], vt = c(5, 4, 5, 5), vm = "Mh", 
                   vh = "Poisson", vtheta = 1.5)
  # The estimators of the demographic parameters obtained with the robust design 
  # are similar to those obtained with the Jolly-Seber model applied to the 
  # between primary period data.
cbind(op.m1$survivals, rd.m4$survivals)
cbind(op.m1$N, rd.m4$N)
cbind(op.m1$birth, rd.m4$birth)
cbind(op.m1$Ntot, rd.m4$Ntot)
}
\keyword{datasets}
