\name{predict.fda}
\alias{predict.fda}
\title{
make predictions from an fda object
}
\usage{
predict.fda(object, x, type, prior, dimension)
}
\arguments{
\item{object}{
An object of class \code{fda}.
}
\item{x}{
new data at which to make predictions. If missing, the training data is used.
}
\item{type}{
kind of predictions: \code{type = class} (default) produces a fitted factor,
\code{type = variates} produces a matrix of discriminant variables,
\code{type = posterior} produces a matrix of posterior probabilities (based on a
gaussian assumption), and \code{type=hierarchical} produces the predicted class in
sequence for models of all dimensions.
}
\item{prior}{
the prior probabability vector for each class; the default is the training 
sample proportions.
}
\item{dimension}{
the dimension of the space to be used, no larger than the dimension component
 of \code{object}.
}}
\value{
An appropriate object depending on \code{type}. \code{object} has a component 
\code{"fit}" which is regression fit produced by the \code{method} argument to \code{fda}.
There should be a \code{predict} method for this object which is invoked. 
This method should itself take as input \code{object} and optionally \code{x}.
}
\seealso{
\code{\link{fda}}, \code{\link{mars}}, \code{\link{bruto}}, \code{\link{polyreg}}, \code{\link{softmax}}, \code{\link{confusion}}
}
% Converted by Sd2Rd version 0.2-a5.
