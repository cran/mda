\name{predict.mda}
\alias{predict.mda}
\title{Classify by Mixture Discriminant Analysis}
\description{
  Classify observations in conjunction with \code{mda}.
}
\usage{
predict.mda(object, x, type, prior, dimension, \dots)
}
\arguments{
  \item{object}{a fitted mda object.}
  \item{x}{new data at which to make predictions.  If missing, the
    training data is used.}
  \item{type}{kind of predictions: \code{type = "class"} (default)
    produces a fitted factor, \code{type = "variates"} produces a matrix
    of discriminant variables (note that the maximal dimension is
    determined by the number of subclasses), \code{type = "posterior"}
    produces a matrix of posterior probabilities (based on a gaussian
    assumption), \code{type = "hierarchical"} produces the predicted
    class in sequence for models of dimensions specified by
    \code{dimension} argument.}
  \item{prior}{the prior probabability vector for each class; the
    default is the training sample proportions.}
  \item{dimension}{the dimension of the space to be used, no larger
    than the dimension component of \code{object}, and in general less
    than the number of subclasses.  \code{dimension} can be a vector for
    use with \code{type = "hierarchical"}.}
}
\value{
  An appropriate object depending on \code{type}.  \code{object} has a
  component \code{fit} which is regression fit produced by the
  \code{method} argument to \code{mda}.  There should be a
  \code{predict} method for this object which is invoked.  This method
  should itself take as input \code{object} and optionally \code{x}.
}
\seealso{
  \code{\link{mda}},
  \code{\link{fda}},
  \code{\link{mars}},
  \code{\link{bruto}},
  \code{\link{polyreg}},
  \code{\link{softmax}},
  \code{\link{confusion}}
}
\examples{
data(glass)
samp <- sample(1:nrow(glass), 100)
glass.train <- glass[samp,]
glass.test <- glass[-samp,]
glass.mda <- mda(Type ~ ., data = glass.train)
predict(glass.mda, glass.test, type = "post") # abbreviations are allowed
confusion(glass.mda, glass.test)
}
\keyword{classif}
% Converted by Sd2Rd version 0.3-3.
