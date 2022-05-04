\name{HIV}
\alias{HIV}
\docType{data}
\title{ Epidemiological capture-recapture Data on HIV}
\description{
  Epidemiological capture-recapture data on HIV from four reporting centers in Rome, Italy.
}
\usage{HIV}
\format{
15 by 5 numeric matrix, with the following columns:
  \describe{
    \item{\code{c1}, \code{c2}, \code{c3}, \code{c4}}{ Capture histories from the four capture occasions}
    \item{\code{freq}}{ Observed frequencies for each capture history}
    }
}
\details{
  The capture histories were obtained by linking the records of the four reporting centers. 

  Each row of this data set represents an observed capture history followed by its frequency. 
}
\source{
Abeni, D.A., Brancato, G. and Perucci, C. A. (1994) Capture-recapture to estimate the size of the population with human immunodeficiency virus type 1 infection.  \emph{Epidemiology}, \bold{5}, 410--414.
}
\references{
Baillargeon, S. and Rivest, L.P. (2007) Rcapture: Loglinear models for capture-recapture in R.  \emph{Journal of Statistical Software}, \bold{19}(5), \doi{10.18637/jss.v019.i05}.
}
\examples{
desc <- descriptive(HIV, dfreq=TRUE)
desc
  # 1774 out of 1896 individuals (94\%) appear on one list only.

plot(desc)
  # The fi plot is linear showing that heterogeneity is not a problem.  
  # Models with a time (or list) effect and possible pairwise dependencies 
  # between lists will be considered.

cp.m1 <- closedpCI.t(HIV, dfreq = TRUE, mX = ~ (c1+c2+c3+c4)^2, 
                     mname = "Mt double interaction")
cp.m1
  # The model fits well. Let's find out which interactions are important.
summary(cp.m1$fit)$coefficients
  # Eliminating the non significant interactions stepwise shows that only 
  # the [1,2] interaction is important.

closedpCI.t(HIV, dfreq = TRUE, mX = ~ . + c1:c2, mname = "Mt interaction 1,2")
}
\keyword{datasets}
