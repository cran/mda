\name{mda.start}
\alias{mda.start}
\title{Initialization for Mixture Discriminant Analysis.}
\description{Provide starting weights for the \code{mda} function which
  performs discriminant analysis by gaussian mixtures.}
\usage{
mda.start(x, g, subclasses = 3, trace.mda.start = FALSE,
          start.method = c("kmeans", "lvq"), tries = 5,
          criterion = c("misclassification", "deviance"), \dots)
}
\arguments{
  \item{x}{The x data, or an mda object.}
  \item{g}{The response vector g.}
  \item{subclasses}{number of subclasses per class, as in \code{mda}.}
  \item{trace.mda.start}{Show results of each iteration.}
  \item{start.method}{Either \code{"kmeans"} or \code{"lvq"}.  The
    latter requires package \bold{class} (from the \bold{VR} package
    bundle.}
  \item{tries}{Number of random starts.}
  \item{criterion}{By default, classification errors on the training
    data.  Posterior deviance is also an option.}
}
\value{
  A list of weight matrices, one for each class.
}
% Converted by Sd2Rd version 0.3-3.
\keyword{classif}
