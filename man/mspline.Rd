\name{mspline}
\alias{mspline}
\title{
Vector Smoothing Spline
}
\description{
Fit a smoothing spline to a matrix of responses, single x.
}
\usage{
mspline(x, y, w, df = 5, lambda, thresh = 1e-04, \dots)
}
\arguments{
  \item{x}{
x variable (numeric vector).
}
  \item{y}{
response matrix.
}
  \item{w}{
optional weight vector, defaults to a vector of ones.
}
  \item{df}{
requested degrees of freedom, as in \code{smooth.spline}.
}
  \item{lambda}{
can provide penalty instead of df.
}
  \item{thresh}{
convergence threshold for df inversion (to lambda).
}
  \item{\dots}{
holdall for other arguments.
}
}
\details{
  This function is based on the ingredients of \code{smooth.spline},
  and allows for simultaneous smoothing of multiple responses
}
\value{
  A list is returned, with  a number of components, only some of which
  are of interest. These are
  \item{lambda}{The value of lambda used (in case df was supplied)}
  \item{df}{The df used (in case lambda was supplied)}
  \item{s}{A matrix like \code{y} of smoothed responses}
  \item{lev}{Self influences (diagonal of smoother matrix)}
}
\author{
Trevor Hastie
}
\seealso{
\code{\link{smooth.spline}}
}
\examples{
x=rnorm(100)
y=matrix(rnorm(100*10),100,10)
fit=mspline(x,y,df=5)
}
