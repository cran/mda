\name{model.matrix.mars}
\alias{model.matrix.mars}
\title{
    produce a model matrix from a `mars' object
}
\usage{
model.matrix.mars(object, x, which, full=F, \dots)
}
\arguments{
    \item{object}{a \code{mars} object}
    \item{x}{optional argument; if supplied, the mars basis functions
	are evaluated at these new observations.}
    \item{which}{which columns should be used. The default is to use the
	columns described by the component \code{selected.terms} on
	\code{object}.}
    \item{full}{if \code{TRUE} the entire set of columns are selected,
	even redundant ones.  This is used for updating a mars fit.}
}
\value{
    a model matrix corresponding to the selected columns
}
\seealso{
    \code{\link{mars}},
    \code{\link{predict.mars}}
}
% Converted by Sd2Rd version 0.3-3.
