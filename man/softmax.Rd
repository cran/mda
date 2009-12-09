\name{softmax}
\alias{softmax}
\title{Find the Maximum in Each Row of a Matrix}
\description{Find the maximum in each row of a matrix.}
\usage{
softmax(x, gap = FALSE)
}
\arguments{
  \item{x}{a numeric matrix.}
  \item{gap}{if \code{TRUE}, the difference between the largest and
    next largest column is returned.}
}
\value{
  A factor with levels the column labels of \code{x} and values the
  columns corresponding to the maximum column. If \code{gap = TRUE} a
  list is returned, the second component of which is the difference
  between the largest and next largest column of \code{x}.
}
\seealso{
  \code{\link{predict.fda}},
  \code{\link{confusion}},
  \code{\link{fda}} 
  \code{\link{mda}} 
}
\examples{
data(iris)
irisfit <- fda(Species ~ ., data = iris)
posteriors <- predict(irisfit, type = "post")
confusion(softmax(posteriors), iris[, "Species"])
}
% Converted by Sd2Rd version 0.3-3.
\keyword{utilities}
