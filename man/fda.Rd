\name{fda}
\alias{fda}
\alias{coef.fda}
\alias{plot.fda}
\alias{print.fda}
\title{Flexible Discriminant Analysis}
\description{
  Flexible discriminant analysis.}
\usage{
fda(formula, data, weights, theta, dimension, eps, method,
    keep.fitted, \dots)
}
\arguments{
  \item{formula}{of the form \code{y~x} it describes the response and
    the predictors.  The formula can be more complicated, such as
    \code{y~log(x)+z} etc (see \code{\link{formula}} for more details).
    The response should be a factor representing the response variable,
    or any vector that can be coerced to such (such as a logical
    variable).}
  \item{data}{data frame containing the variables in the formula
    (optional).}
  \item{weights}{an optional vector of observation weights.}
  \item{theta}{an optional matrix of class scores, typically with less
    than \code{J-1} columns.}
  \item{dimension}{The dimension of the solution, no greater than
    \code{J-1}, where \code{J} is the number classes.  Default is
    \code{J-1}.}
  \item{eps}{a threshold for small singular values for excluding
    discriminant variables; default is \code{.Machine$double.eps}.}
  \item{method}{regression method used in optimal scaling.  Default is
    linear regression via the function \code{polyreg}, resulting in
    linear discriminant analysis.  Other possibilities are \code{mars}
    and \code{bruto}.  For Penalized Discriminant analysis
    \code{gen.ridge} is appropriate.}
  \item{keep.fitted}{a logical variable, which determines whether the
    (sometimes large) component \code{"fitted.values"} of the \code{fit}
    component of the returned fda object should be kept.  The default is
    \code{TRUE} if \code{n * dimension < 1000}.}
  \item{\dots}{additional arguments to \code{method}.}
}
\value{
  an object of class \code{"fda"}.  Use \code{predict} to extract
  discriminant variables, posterior probabilities or predicted class
  memberships.  Other extractor functions are \code{coef},
  \code{confusion} and \code{plot}.

  The object has the following components:
  \item{percent.explained}{the percent between-group variance explained
    by each dimension (relative to the total explained.)}
  \item{values}{optimal scaling regresssion sum-of-squares for each
    dimension (see reference).  The usual discriminant analysis
    eigenvalues are given by \code{values / (1-values)}, which are used
    to define \code{percent.explained}.}
  \item{means}{class means in the discriminant space.  These are also
    scaled versions of the final theta's or class scores, and can be
    used in a subsequent call to \code{fda} (this only makes sense if
    some columns of theta are omitted---see the references).} 
  \item{theta.mod}{(internal) a class scoring matrix which allows
    \code{predict} to work properly.}
  \item{dimension}{dimension of discriminant space.}
  \item{prior}{class proprotions for the training data.}
  \item{fit}{fit object returned by \code{method}.}
  \item{call}{the call that created this object (allowing it to be
    \code{update}-able)}
  \item{confusion}{confusion matrix when classifying the training data.}

  The \code{method} functions are required to take arguments \code{x}
  and \code{y} where both can be matrices, and should produce a matrix
  of \code{fitted.values} the same size as \code{y}.  They can take
  additional arguments \code{weights} and should all have a \code{\dots}
  for safety sake.  Any arguments to \code{method} can be passed on via
  the \code{\dots} argument of \code{fda}.  The default method
  \code{\link{polyreg}} has a \code{degree} argument which allows
  polynomial regression of the required total degree.  See the
  documentation for \code{\link{predict.fda}} for further requirements
  of \code{method}.
}
\note{
  This software it is not well-tested, we would like to hear of any
  bugs.}
\author{
  Trevor Hastie and Robert Tibshirani
}
\seealso{
  \code{\link{predict.fda}},
  \code{\link{mars}},
  \code{\link{bruto}},
  \code{\link{polyreg}},
  \code{\link{softmax}},
  \code{\link{confusion}},
  %% \code{\link{coefficients.fda}}
}
\references{
  ``Flexible Disriminant Analysis by Optimal Scoring''  by Hastie,
  Tibshirani and Buja, 1994, JASA, 1255-1270.

  ``Penalized Discriminant Analysis'' by Hastie, Buja and Tibshirani,
  Annals of Statistics, 1995 (in press).
}
\examples{
data(iris)
irisfit <- fda(Species ~ ., data = iris)
irisfit
## fda(formula = Species ~ ., data = iris)
##
## Dimension: 2 
##
## Percent Between-Group Variance Explained:
##     v1     v2 
##  99.12 100.00 
##
## Degrees of Freedom (per dimension): 5 
##
## Training Misclassification Error: 0.02 ( N = 150 )

confusion(irisfit, iris)
##            Setosa Versicolor Virginica 
##     Setosa     50          0         0
## Versicolor      0         48         1
##  Virginica      0          2        49
## attr(, "error"):
## [1] 0.02

plot(irisfit)

coef(irisfit)
##           [,1]        [,2]
## [1,] -2.126479 -6.72910343
## [2,] -0.837798  0.02434685
## [3,] -1.550052  2.18649663
## [4,]  2.223560 -0.94138258
## [5,]  2.838994  2.86801283

marsfit <- fda(Species ~ ., data = iris, method = mars)
marsfit2 <- update(marsfit, degree = 2)
marsfit3 <- update(marsfit, theta = marsfit$means[, 1:2]) 
## this refits the model, using the fitted means (scaled theta's)
## from marsfit to start the iterations
}
\keyword{classif}
% Converted by Sd2Rd version 0.3-3.
