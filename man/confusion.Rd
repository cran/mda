\name{confusion}
\alias{confusion}
\title{
    compute the confusion matrix between two factors, or for an fda or
    mda object
}
\usage{
confusion(predict, true)
}
\arguments{
    \item{predict}{the predicted factor, or an fda or mda model}
    \item{true}{the true factor, or else (in the case where predict is a
	model object) a data frame (list) containing the test data.}
}
\value{
    essential table(predict, true), but with some useful attribute(s)
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
