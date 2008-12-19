\name{model.matrix.mars}
\alias{model.matrix.mars}
\title{Produce a Design Matrix from a `mars' Object}
\description{Produce a design matrix from a `mars' object.}
\usage{
\method{model.matrix}{mars}(object, x, which, full = FALSE, \dots)
}
\arguments{
  \item{object}{a mars object.}
  \item{x}{optional argument; if supplied, the mars basis functions are
    evaluated at these new observations.}
  \item{which}{which columns should be used. The default is to use the
    columns described by the component \code{selected.terms} on
    \code{object}.}
  \item{full}{if \code{TRUE} the entire set of columns are selected,
    even redundant ones.  This is used for updating a mars fit.}
  \item{\dots}{further arguments to be passed from or to methods.}
}
\value{
  A model matrix corresponding to the selected columns.
}
\seealso{
  \code{\link{mars}},
  \code{\link{predict.mars}}
}
\keyword{models}
