\name{predict.fda}
\alias{predict.fda}
\title{
    make predictions from an fda object
}
\usage{
predict.fda(object, x, type, prior, dimension)
}
\arguments{
    \item{object}{An object of class \code{fda}.}
    \item{x}{New data at which to make predictions. If missing, the
	training data is used.}
    \item{type}{kind of predictions: \code{type = class} (default)
	produces a fitted factor, \code{type = variates} produces a
	matrix of discriminant variables, \code{type = posterior}
	produces a matrix of posterior probabilities (based on a
	gaussian assumption), and \code{type=hierarchical} produces the
	predicted class in sequence for models of all dimensions.}
    \item{prior}{the prior probabability vector for each class; the
	default is the training sample proportions.}
    \item{dimension}{the dimension of the space to be used, no larger
	than the dimension component  of \code{object}.}
}
\value{
    An appropriate object depending on \code{type}. \code{object} has a
    component \code{"fit}"} which is regression fit produced by the
    \code{method} argument to \code{fda}.  There should be a
    \code{predict} method for this object which is invoked.  This method
    should itself take as input \code{object} and optionally \code{x}.
}
\seealso{
    \code{\link{fda}},
    \code{\link{mars}},
    \code{\link{bruto}},
    \code{\link{polyreg}},
    \code{\link{softmax}},
    \code{\link{confusion}}
}
\examples{
data(iris)
irisfit <- fda(Species ~ ., data = iris)
irisfit
## Call:
## fda(x = iris$x, g = iris$g)
## 
## Dimension: 2 
##
## Percent Between-Group Variance Explained:
##     v1  v2 
##  99.12 100
confusion(predict(irisfit, iris), iris$Species)
##            Setosa Versicolor Virginica
##     Setosa     50          0         0
## Versicolor      0         48         1
##  Virginica      0          2        49
## attr(, "error"):
## [1] 0.02
}
% Converted by Sd2Rd version 0.3-3.
