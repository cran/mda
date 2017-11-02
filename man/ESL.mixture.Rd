\name{ESL.mixture}
\alias{ESL.mixture}
\title{Mixture example from "Elements of Statistical Learning"}
\description{
  A list with training data and other details for the mixture example
}
\usage{data(ESL.mixture)}
\format{
  This list contains the following elements:
  \describe{
    \item{x}{a 200x2 matrix of predictors.}
    \item{y}{a 200 vector of y values taking values 0 or 1.}
    \item{xnew}{a 6831x2 matrix of prediction points, on a 69x99 grid.}
    \item{prob}{a vector of 6831 probabilities - the true probabilities
      of a 1 at each point in \code{xnew}.}
    \item{marginal}{the marginal distribution of the predictors t each
      point in \code{xnew}.}
    \item{px1}{grid values for first coordinate in \code{xnew}.}
    \item{px2}{grid values for second coordinate in \code{xnew}.}
    \item{means}{a 20 x 2 matrix of means used in the generation of
      these data.}
  }
}
\source{
  "Elements of Statistical Learning (second edition)", Hastie, T., Tibshirani, R. and
  Friedman, J. (2009), Springer, New York.
  \url{http://www.stanford.edu/~hastie/ElemStatLearn}
}
\keyword{datasets}
