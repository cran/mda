\name{predict.bruto}
\alias{predict.bruto}
\title{Predict method for BRUTO Objects}
\description{
  Predicted values based on `bruto' additive spline models which are fit
  by adaptive backfitting.
}
\usage{
\method{predict}{bruto}(object, x, type=c("fitted", "terms"), \dots)
}
\arguments{
  \item{object}{a fitted bruto object}
  \item{x}{values at which predictions are to be made.}
  \item{type}{if type is \code{"fitted"}, the fitted values are
    returned.  If type is \code{"terms"}, a list of fitted terms is
    returned, each with an \code{x} and \code{y} component.  These can
    be used to show the fitted functions.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\value{
  Either a fit matrix or a list of fitted terms.
}
\seealso{
  \code{\link{bruto}},
  \code{\link{predict}}
}
\examples{
data(trees)
fit1 <- bruto(trees[,-3], trees[3])
fitted.terms <- predict(fit1, as.matrix(trees[,-3]), type = "terms")
par(mfrow=c(1,2), pty="s")
for(tt in fitted.terms) plot(tt, type="l")
}
% Converted by Sd2Rd version 0.3-3.
\keyword{smooth}
