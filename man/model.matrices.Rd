\name{model.matrices}
\alias{Zdelta}
\alias{Xomega.t}
\alias{Xomega.0}
\alias{Xclosedp}
\title{ Construction of Matrices Useful for Loglinear Models in Capture-Recature Experiments (internal use) }
\description{
\code{Zdelta} produces the Zdelta matrix, i.e. the first component of the design matrix in a loglinear model for open population and robust design models.

\code{Xomega.t} and \code{Xomega.0} produce Xomega matrices, i.e. the second component of the design matrix in a loglinear model for robust designs.
\code{Xomega.t} is used by \code{robustd.t}, and \code{Xomega.0} is used by \code{robustd.0}.

\code{Xclosedp} produces the design matrix for a closed population model. This matrix is used by \code{robustd.t} and \code{robustd.0} to compute the capture probabilities.
}
\usage{
Zdelta(Xdelta)
Xomega.t(vt, vm, vh, va, rd.call)
Xomega.0(vt, vm, vh, va, rd.call)
Xclosedp(t, m, h, a)
}
\arguments{
  \item{Xdelta}{ The Xdelta matrix, i.e. the second component of the the design matrix in loglinear models for open populations or for the robust design.}
  \item{vt}{ A vector containing the numbers of capture occasions for each primary sampling period. The length of this vector equals the number of primary sampling periods (noted \eqn{I}). }
  \item{vm}{ A vector indicating the closed population model for each primary period. The elememts of \code{vm} can be "none"=no model, "M0"=M0 model, "Mt"=Mt model, "Mh"=Mh model or "Mth"=Mth model for \code{Xomega.t}. For \code{Xomega.0}, models with time effect are not allowed. So \code{vm} can than be "none"=no model, "M0"=M0 model, or "Mh"=Mh model. }
  \item{vh}{ A list indicating, for each primary period with a heterogeneity model, the form of the columns for heterogeneity in the design matrix. The elements of \code{vht} can be "Chao", "Poisson", "Darroch" or any numerical \code{R} function created beforehand by the user. "Chao" represents Chao's model, "Poisson" represents the function \eqn{f(k)=a^k-1}, where \eqn{k} is the number of captures, and "Darroch" represents the function \eqn{f(k)=k^2/2}. If an \code{R} function is given, it is the implemantation of any mathematical function \eqn{f(k)}. It has only one argument. For the Poisson model, the exponent's base \code{a} is specified in the argument \code{va}. }
  \item{va}{ A vector indicating, for each primary period with a Poisson model, the value of the exponent's base in \eqn{f(k)=a^k-1}. }
  \item{rd.call}{ The call of the \code{robust} function in which the \code{Xomega} function is called.}
  \item{t}{ The number of capture occasions. }
  \item{m}{ A character string identifying the capture-recapture model : "none"=no model, "M0"=M0 model, "Mt"=Mt model, "Mh"=Mh model or "Mth"=Mth model.}
  \item{h}{ A character string ("Chao", "Poisson" or "Darroch") or a numerical \code{R} function specifying the form of the column for heterogeneity in the design matrix. "Chao" represents Chao's model, "Poisson" represents the function \eqn{f(k)=a^k-1}, where \eqn{k} is the number of captures, and "Darroch" represents the function \eqn{f(k)=k^2/2}. If an \code{R} function is given, it is the implemantation of any mathematical function \eqn{f(k)}. It has only one argument.}
  \item{a}{ The value of the exponent's base for a Poisson model. }
}
\details{
These functions are internal. They are not intended to be called by the user. No argument validation is done within these functions.
}
\value{
\code{Zdelta} returns a matrix with \eqn{2(I-1)} columns and as many rows as the Xdelta matrix.\cr
\code{Xomega.t} and \code{Xomega.0} return a list containing the first four elements of the following enumeration.\cr \code{Xclosedp} returns a list containing the last two elements of the following enumeration.
\item{mat}{ The Xomega matrix containing the columns of the design matrices for the closed population models in every
primary period of a robust design.}
\item{nbparam}{ A vector of length \eqn{I} containing the numbers of parameters for the closed population models
within each primary period.}
\item{models}{ A vector of length \eqn{I} identifying the closed population models for each primary period.}
\item{paramnames}{ A vector of character strings identifying the parameters.}
\item{mat}{ The design matrix for the closed population model.}
\item{nbparam}{ The number of parameters in the closed population model.}
}
\references{ 
Rivest, L.P. and Daigle, G. (2004) Loglinear models for the robust design in mark-recapture experiments. \emph{Biometrics}, \bold{60}, 100-107.
}
\author{Sophie Baillargeon \email{sbaillar@mat.ulaval.ca} and Louis-Paul Rivest \email{lpr@mat.ulaval.ca} }
\note{
The \code{Zdelta} function is called by \code{\link{openp}}, \code{\link{robustd.t}}, \code{\link{robustd.0}}.

The \code{Xomega.t} and \code{Xclosedp} functions are called by \code{\link{robustd.t}}.

The \code{Xomega.0} and \code{Xclosedp} function is called by \code{\link{robustd.0}}.
}
\keyword{ manip }
\keyword{ internal }
