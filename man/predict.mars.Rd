\name{predict.mars}
\alias{predict.mars}
\title{Predict method for MARS Objects}
\description{
  Predicted values based on `mars' multivariate adaptive regression
  spline models.
}
\usage{
\method{predict}{mars}(object, x, \dots)
}
\arguments{
  \item{object}{an object of class \code{"mars"}.}
  \item{x}{values at which predictions are to be made.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\value{
  the fitted values.
}
\seealso{
  \code{\link{mars}},
  \code{\link{predict}},
  \code{\link{model.matrix.mars}}
}
%%\examples{
%%## see help for mars()
%%}
% Converted by Sd2Rd version 0.3-3.
\keyword{smooth}
