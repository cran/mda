\name{gen.ridge}
\alias{gen.ridge}
\alias{predict.gen.ridge}
\title{Penalized Regression}
\description{
  Perform a penalized regression, as used in penalized discriminant
  analysis.
}
\usage{
gen.ridge(x, y, weights, lambda=1, omega, df, \dots)
}
\arguments{
  \item{x, y, weights}{the x and y matrix and possibly a weight vector.}
  \item{lambda}{the shrinkage penalty coefficient.}
  \item{omega}{a penalty object; omega is the eigendecomposition of
    the penalty matrix, and need not have full rank.  By default,
    standard ridge is used.}
  \item{df}{an alternative way to prescribe lambda, using the notion
    of equivalent degrees of freedom.}
  \item{\dots}{currently not used.}
}
\value{
  A generalized ridge regression, where the coefficients are penalized
  according to omega.  See the function definition for further details.
  No functions are provided for producing one dimensional penalty
  objects (omega).
  %%FIXME
  %% \code{\link{laplacian}()} creates a two dimensional penalty
  %% object, suitable for (small) images.
}
%% \seealso{
%%   \code{\link{laplacian}}
%% }
% Converted by Sd2Rd version 0.3-3.
\keyword{regression}
