\name{laplacian}
\Rdversion{1.1}
\alias{laplacian}
\title{create penalty object for two-dimensional smoothing.
}
\description{
  Creates a penalty matrix  for use by \code{gen.ridge} for
  two-dimensional smoothing.
}
\usage{
laplacian(size, compose)
laplacian(size = 16, compose = FALSE)
}
\arguments{
  \item{size}{dimension of the image is \code{size x size}; default is 16.}
  \item{compose}{default is \code{compose=FALSE}, which means the
    penalty is returned as an eigen-decomposition. If \code{compose=TRUE}, a penalty matrix is returned.}
}
\details{
Formulas are used to construct a laplacian for smoothing a square image.
}
\value{
  If \code{compose=FALSE}, an eigen-decomposition object is
  returned. The \code{vectors} component is a \code{size^2 x size^2}
  orthogonal matrix, and the \code{$values} component is a \code{size^2}
  vector of non-negative eigen-values. If \code{compose=TRUE}, these are
  multiplied together to form a single matrix.
}
\references{
Here we follow very closely the material on page 635 in JASA 1991 of O'Sullivan's article on discretized Laplacian Smoothing
}
\author{
Trevor Hastie <hastie@stanford.edu
}
\seealso{\code{\link{gen.ridge}},\code{\link{fda}}
}
\keyword{regression}

