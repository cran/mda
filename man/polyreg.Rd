\name{polyreg}
\alias{polyreg}
\alias{predict.polyreg}
\title{Polynomial Regression}
\description{
  Simple minded polynomial regression.
}
\usage{
polyreg(x, y, w, degree = 1, monomial = FALSE, \dots)
}
\arguments{
  \item{x}{predictor matrix.}
  \item{y}{response matrix.}
  \item{w}{optional (positive) weights.}
  \item{degree}{total degree of polynomial basis (default is 1).}
  \item{monomial}{If \code{TRUE} a monomial basis is used (no cross
    terms).  Default is \code{FALSE}.}
  \item{\dots}{currently not used.}
}
\value{
  A polynomial regression fit, containing the essential ingredients for
  its \code{predict} method.
}
%%\seealso{
%%  \code{\link{predict.polyreg}}
%%}
% Converted by Sd2Rd version 0.3-3.
\keyword{regression}
