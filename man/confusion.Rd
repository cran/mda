\name{confusion}
\alias{confusion}
\alias{confusion.default}
\alias{confusion.list}
\title{Confusion Matrices}
\description{
  Compute the confusion matrix between two factors, or for an fda or
  mda object.
}
\usage{
\method{confusion}{default}(predict, true, \dots)
\method{confusion}{fda}(object, data, \dots)
}
\arguments{
  \item{predict}{the predicted factor.}
  \item{true}{the true factor.}
  \item{object}{an fda or mda model object.}
  \item{data}{a data frame (list) containing the test data.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\details{
  This is a generic function.
}
\value{
  For the default method essentially \code{table(predict, true)}, but
  with some useful attribute(s).
}
\seealso{
  \code{\link{fda}},
  \code{\link{predict.fda}}
}
\examples{
data(iris)
irisfit <- fda(Species ~ ., data = iris)
confusion(predict(irisfit, iris), iris$Species)
##            Setosa Versicolor Virginica 
##     Setosa     50          0         0
## Versicolor      0         48         1
##  Virginica      0          2        49
## attr(, "error"):
## [1] 0.02
}
% Converted by Sd2Rd version 0.3-3.
\keyword{category}
