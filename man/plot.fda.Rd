\name{plot.fda}
\alias{plot.fda}
\title{Plot for Flexible Discriminant Analysis}
\description{Plot in discriminant coordinates a \code{fda} or (by inheritance) a \code{mda} object.
}
\usage{
\method{plot}{fda}(x, data, coords, group, colors, pch, mcolors, mpch, pcex, mcex, \dots)
}
\arguments{
  \item{x}{an object of class \code{"fda"}.}
  \item{data}{the data to plot in the discriminant coordinates. If
    \code{group="true"}, then data should be a data frame with the same variables that were used in the fit. If \code{group="predicted"},
    \code{data} need not contain the response variable, and can in fact
    be the correctly-sized \code{"x"} matrix.}
  \item{coords}{vector of coordinates to plot, with default
    \code{coords="c(1,2)"}. All pairs of plots are produced.}
  \item{group}{if \code{group="true"} (the default), each point is color
    and symbol coded according to the response in \code{data}. If
    \code{group="predicted"}, the class of each point is predicted from
    the model, and used instead.}
  \item{colors}{a vector of colors to be used in the plotting.}
  \item{pch}{a vector of plotting characters.}
  \item{mcolors}{a vector of colors for the class centroids; default is \code{colors}.}
  \item{mpch}{a vector of plotting characters for the centroids.}
  \item{pcex}{character expansion factor for the points; defualt is \code{pcex="0.5"}.}
  \item{mcex}{character expansion factor for the centroids; defualt is \code{pcex="2.5"}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\seealso{
   \code{\link{fda}},
   \code{\link{mda}},
   \code{\link{predict.fda}}
 }
\examples{
data(iris)
irisfit <- fda(Species ~ ., data = iris)
plot(irisfit)
data(ESL.mixture)
## Not a data frame
mixture.train=ESL.mixture[c("x","y")] 
mixfit=mda(y~x, data=mixture.train)
plot(mixfit, mixture.train)
plot(mixfit, data=ESL.mixture$xnew, group="pred")
}
\keyword{classif}
