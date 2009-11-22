\name{mda}
\alias{mda}
\alias{print.mda}
\title{Mixture Discriminant Analysis}
\description{
  Mixture discriminant analysis.
}
\usage{
mda(formula, data, subclasses, sub.df, tot.df, dimension, eps,
    iter, weights, method, keep.fitted, trace, \dots)
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
  \item{subclasses}{Number of subclasses per class, default is 3.  Can be
    a vector with a number for each class.}
  \item{sub.df}{If subclass centroid shrinking is performed, what is the
    effective degrees of freedom of the centroids per class.  Can be a
    scalar, in which case the same number is used for each class, else a
    vector.}
  \item{tot.df}{The total df for all the centroids can be specified
    rather than separately per class.}
  \item{dimension}{The dimension of the reduced model.  If we know our
    final model will be confined to a discriminant subspace (of the
    subclass centroids), we can specify this in advance and have the EM
    algorithm operate in this subspace.}
  \item{eps}{A numerical threshold for automatically truncating the
    dimension.}
  \item{iter}{A limit on the total number of iterations,  default is 5.}
  \item{weights}{\emph{NOT} observation weights!  This is a special
    weight structure, which for each class assigns a weight (prior
    probability) to each of the observations in that class of belonging
    to one of the subclasses.  The default is provided by a call to
    \code{mda.start(x, g, subclasses, trace, \dots)} (by this time
    \code{x} and \code{g} are known).  See the help for
    \code{\link{mda.start}}.  Arguments for \code{mda.start} can be
    provided via the \code{\dots} argument to mda, and the
    \code{weights} argument need never be accessed.  A previously fit
    mda object can be supplied, in which case the final subclass
    \code{responsibility} weights are used for \code{weights}.  This 
    allows the iterations from a previous fit to be continued.}
  \item{method}{regression method used in optimal scaling.  Default is
    linear regression via the function \code{polyreg}, resulting in the
    usual mixture model.  Other possibilities are \code{mars} and 
    \code{bruto}.  For penalized mixture discriminant models
    \code{gen.ridge} is appropriate.}
  \item{keep.fitted}{a logical variable, which determines whether the
    (sometimes large) component \code{"fitted.values"} of the \code{fit}
    component of the returned \code{mda} object should be kept.  The
    default is \code{TRUE} if \code{n * dimension < 1000}.} 
  \item{trace}{if \code{TRUE}, iteration information is printed.  Note
    that the deviance reported is for the posterior class likelihood,
    and not the full likelihood, which is used to drive the EM algorithm
    under \code{mda}.  In general the latter is not available.}
  \item{\dots}{additional arguments to \code{mda.start} and to
    \code{method}.}
}
\value{
  An object of class \code{c("mda", "fda")}.  The most useful extractor
  is \code{predict}, which can make many types of predictions from this
  object.  It can also be plotted, and any functions useful for fda
  objects will work here too, such as \code{confusion} and \code{coef}.

  The object has the following components:
  \item{percent.explained}{the percent between-group variance explained
    by each dimension (relative to the total explained.)}
  \item{values}{optimal scaling regression sum-of-squares for each
    dimension (see reference).}
  \item{means}{subclass means in the discriminant space.  These are also
    scaled versions of the final theta's or class scores, and can be
    used in a subsequent call to \code{mda} (this only makes sense if
    some columns of theta are omitted---see the references)}
  \item{theta.mod}{(internal) a class scoring matrix which allows
    \code{predict} to work properly.}
  \item{dimension}{dimension of discriminant space.}
  \item{sub.prior}{subclass membership priors, computed in the fit.  No
    effort is currently spent in trying to keep these above a threshold.}
  \item{prior}{class proportions for the training data.}
  \item{fit}{fit object returned by \code{method}.}
  \item{call}{the call that created this object (allowing it to be
    \code{update}-able).}
  \item{confusion}{confusion matrix when classifying the training data.}
  \item{weights}{These are the subclass membership probabilities for
    each member of the training set; see the weights argument.}
  \item{assign.theta}{a pointer list which identifies which elements of
    certain lists belong to individual classes.}
  \item{deviance}{The multinomial log-likelihood of the fit.  Even though
    the full log-likelihood drives the iterations, we cannot in general
    compute it because of the flexibility of the \code{method} used.
    The deviance can increase with the iterations, but generally does not.}

  The \code{method} functions are required to take arguments \code{x}
  and \code{y} where both can be matrices, and should produce a matrix
  of \code{fitted.values} the same size as \code{y}.  They can take
  additional arguments \code{weights} and should all have a \code{\dots}
  for safety sake.  Any arguments to method() can be passed on via the
  \code{\dots} argument of \code{mda}.  The default method
  \code{polyreg} has a \code{degree} argument which allows polynomial
  regression of the required total degree.  See the documentation for
  \code{\link{predict.fda}} for further requirements of \code{method}.

  The function \code{mda.start} creates the starting weights; it takes
  additional arguments which can be passed in via the \code{\dots}
  argument to \code{mda}.  See the documentation for \code{mda.start}.
}
\note{
  This software it is not well-tested, we would like to hear of any
  bugs.
}
\author{
  Trevor Hastie and Robert Tibshirani
}
\seealso{
  \code{\link{predict.mda}},
  \code{\link{mars}},
  \code{\link{bruto}},
  \code{\link{polyreg}},
  \code{\link{gen.ridge}},
  \code{\link{softmax}},
  \code{\link{confusion}}
  %%\code{\link{coef.fda}},
  %%\code{\link{plot.fda}}
}
\references{
  ``Flexible Disriminant Analysis by Optimal Scoring'' by Hastie,
  Tibshirani and Buja, 1994, JASA, 1255-1270.

  ``Penalized Discriminant Analysis'' by Hastie, Buja and Tibshirani,
  Annals of Statistics, 1995 (in press).
    
  ``Discriminant Analysis by Gaussian Mixtures'' by Hastie and
  Tibshirani, 1994, JRSS-B (in press).
}
\examples{
data(iris)
irisfit <- mda(Species ~ ., data = iris)
irisfit
## Call:
## mda(formula = Species ~ ., data = iris)
##
## Dimension: 4
##
## Percent Between-Group Variance Explained:
##     v1     v2     v3     v4
##  96.02  98.55  99.90 100.00
##
## Degrees of Freedom (per dimension): 5
##
## Training Misclassification Error: 0.02 ( N = 150 )
##
## Deviance: 15.102

data(glass)
# random sample of size 100
samp <- c(1, 3, 4, 11, 12, 13, 14, 16, 17, 18, 19, 20, 27, 28, 31,
          38, 42, 46, 47, 48, 49, 52, 53, 54, 55, 57, 62, 63, 64, 65,
          67, 68, 69, 70, 72, 73, 78, 79, 83, 84, 85, 87, 91, 92, 94,
          99, 100, 106, 107, 108, 111, 112, 113, 115, 118, 121, 123,
          124, 125, 126, 129, 131, 133, 136, 139, 142, 143, 145, 147,
          152, 153, 156, 159, 160, 161, 164, 165, 166, 168, 169, 171,
          172, 173, 174, 175, 177, 178, 181, 182, 185, 188, 189, 192,
          195, 197, 203, 205, 211, 212, 214) 
glass.train <- glass[samp,]
glass.test <- glass[-samp,]
glass.mda <- mda(Type ~ ., data = glass.train)
predict(glass.mda, glass.test, type="post") # abbreviations are allowed
confusion(glass.mda,glass.test)
}
\keyword{classif}
% Converted by Sd2Rd version 0.3-3.
