\name{predict.mars}
\alias{predict.mars}
\title{Predict method for MARS Objects}
\description{
  Predicted values based on `mars' multivariate adaptive regression
  spline models.
}
\usage{
predict.mars(object, x)
}
\arguments{
  \item{object}{an object of class \code{"mars"}.}
  \item{x}{values at which predictions are to be made.}
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
