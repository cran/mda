\name{mda.start}
\alias{mda.start}
\title{
    Provide starting weights for the mda function
}
\usage{
mda.start(x, g, subclasses = 3, trace.mda.start = FALSE,
          start.method = c("kmeans", "lvq"), tries = 5,
          criterion = c("misclassification", "deviance"), \dots)
}
\arguments{
    \item{x}{The x data, or an mda object}
    \item{g}{The response vector g}
    \item{subclasses}{number of subclasses per class, as in mda.}
    \item{trace.mda.start}{Show results of each iteration}
    \item{start.method}{Either "kmeans" or "lvq". The latter requires
	the statlib contribution "classif".}
    \item{tries}{Number of random starts}
    \item{criterion}{By default, classification errors on the training
	data. Posterior deviance is also an option.}
}
\value{
    A list of weight matrices, one for each class
}
% Converted by Sd2Rd version 0.3-3.
