\name{bruto}
\alias{bruto}
\title{
  Fit an Additive Spline Model by Adaptive Backfitting
}
\description{
  Fit an additive spline model by adaptive backfitting.
}
\usage{
bruto(x, y, w, wp, dfmax, cost, maxit.select, maxit.backfit, 
      thresh = 0.0001, trace.bruto = FALSE, start.linear = TRUE,
      fit.object, \dots)
}
\arguments{
  \item{x}{a matrix of numeric predictors (does not include the column
    of 1s).}
  \item{y}{a vector or matrix of responses.}
  \item{w}{optional observation weight vector.}
  \item{wp}{optional weight vector for each column of \code{y}; the RSS
    and GCV criteria use a weighted sum of squared residuals.}
  \item{dfmax}{a vector of maximum df (degrees of freedom) for each
    term.}
  \item{cost}{cost per degree of freedom; default is 2.}
  \item{maxit.select}{maximum number of iterations during the selection
    stage.}
  \item{maxit.backfit}{maximum number of iterations for the final
    backfit stage (with fixed lambda).} 
  \item{thresh}{convergence threshold (default is 0.0001); iterations
    cease when the relative change in GCV is below this threshold.} 
  \item{trace.bruto}{logical flag.  If \code{TRUE} (default) a progress
    report is printed during the fitting.}
  \item{start.linear}{logical flag.  If \code{TRUE} (default), the model
    starts with the linear fit.}
  \item{fit.object}{This the object returned by \code{bruto()}; if
    supplied, the same model is fit to the presumably new \code{y}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\value{
  A multiresponse additive model fit object of class \code{"bruto"} is
  returned.  The model is fit by adaptive backfitting using smoothing
  splines.  If there are \code{np} columns in \code{y}, then \code{np}
  additive models are fit, but the same amount of smoothing (df) is 
  used for  each term.  The procedure chooses between \code{df = 0}
  (term omitted), \code{df = 1} (term linear) or \code{df > 0} (term
  fitted by smoothing spline).  The model selection is based on an
  approximation to the  GCV criterion, which is used at each step of the
  backfitting procedure.  Once the selection process stops, the model is
  backfit using the chosen amount of smoothing.

  A bruto object has the following components of interest:
  \item{lambda}{a vector of chosen smoothing parameters, one for each
    column of \code{x}.}
  \item{df}{the df chosen for each column of \code{x}.}
  \item{type}{a factor with levels \code{"excluded"}, \code{"linear"} or
    \code{"smooth"}, indicating the status of each column of \code{x}.}
  \item{gcv.select gcv.backfit df.select}{The sequence of gcv values and
    df selected during the execution of the function.}
  \item{nit}{the number of iterations used.}
  \item{fitted.values}{a matrix of fitted values.}
  \item{residuals}{a matrix of residuals.}
  \item{call}{the call that produced this object.}
}
\references{
  Trevor Hastie and Rob Tibshirani,
  \emph{Generalized Additive Models},
  Chapman and Hall, 1990 (page 262).

  Trevor Hastie, Rob Tibshirani and Andreas Buja
  ``Flexible Discriminant Analysis by Optimal Scoring''
  JASA 1994, 89, 1255-1270.
}
\seealso{
  \code{\link{predict.bruto}}
}
\examples{
data(trees)
fit1 <- bruto(trees[,-3], trees[3])
fit1$type
fit1$df
## examine the fitted functions
par(mfrow=c(1,2), pty="s")
Xp <- matrix(sapply(trees[1:2], mean), nrow(trees), 2, byrow=TRUE)
for(i in 1:2) {
  xr <- sapply(trees, range)
  Xp1 <- Xp; Xp1[,i] <- seq(xr[1,i], xr[2,i], len=nrow(trees))
  Xf <- predict(fit1, Xp1)
  plot(Xp1[ ,i], Xf, xlab=names(trees)[i], ylab="", type="l")
}
}
% Converted by Sd2Rd version 0.3-3.
\keyword{smooth}
