\name{mars}
\alias{mars}
\title{Multivariate Additive Regression Splines}
\description{
  Multivariate additive regression splines.
}
\usage{
mars(x, y, w, wp, degree, nk, penalty, thresh, prune, trace.mars,
     forward.step, prevfit, \dots)
}
\arguments{
  \item{x}{a matrix containing the independent variables.}
  \item{y}{a vector containing the response variable, or in the case of
    multiple responses, a matrix whose columns are the response values
    for each variable.}
  \item{w}{an optional vector of observation weights.}
  \item{wp}{an optional vector of response weights.}
  \item{degree}{an optional integer specifying maximum interaction
    degree (default is 1).}
  \item{nk}{an optional integer specifying the maximum number of model
    terms.}
  \item{penalty}{an optional value specifying the cost per degree of
    freedom charge (default is 2).}
  \item{thresh}{an optional value specifying forward stepwise stopping
    threshold (default is 0.001).}
  \item{prune}{an optional logical value specifying whether the model
    should be pruned in a backward stepwise fashion (default is
    \code{TRUE}).}
  \item{trace.mars}{an optional logical value specifying whether info
    should be printed along the way (default is \code{FALSE}).}
  \item{forward.step}{an optional logical value specifying whether
    forward stepwise process should be carried out (default is
    \code{TRUE}).}
  \item{prevfit}{optional data structure from previous fit.  To see the
    effect of changing the penalty paramater, one can use prevfit with
    \code{forward.step = FALSE}.}
  \item{\dots}{further arguments to be passed to or from methods.}
}
\value{
  An object of class \code{"mars"}, which is a list with the following
  components:

  \item{call}{call used to \code{mars}.}
  \item{all.terms}{term numbers in full model.  \code{1} is the constant
    term.  Remaining terms are in pairs (\code{2 3}, \code{4 5}, and so
    on).  \code{all.terms} indicates nonsingular set of terms.}
  \item{selected.terms}{term numbers in selected model.}
  \item{penalty}{the input penalty value.}
  \item{degree}{the input degree value.}
  \item{thresh}{the input threshold value.}
  \item{gcv}{gcv of chosen model.}
  \item{factor}{matrix with \eqn{ij}-th element equal to 1 if term
    \eqn{i} has a factor of the form \eqn{x_j > c}, equal to \eqn{-1} if
    term \eqn{i} has a factor of the form \eqn{x_j \le c}, and to 0 if
    \eqn{xj} is not in term \eqn{i}.}
  \item{cuts}{matrix with \eqn{ij}-th element equal to the cut point
    \eqn{c} for variable \eqn{j} in term \eqn{i}.}
  \item{residuals}{residuals from fit.}
  \item{fitted}{fitted values from fit.}
  \item{lenb}{length of full model.}
  \item{coefficients}{least squares coefficients for final model.}
  \item{x}{a matrix of basis functions obtained from the input x
    matrix.}
}
\note{
  This function was coded from scratch, and did not use any of
  Friedman's mars code.  It gives quite similar results to  Friedman's
  program in our tests, but not exactly the same results.  We have not 
  implemented Friedman's anova decomposition nor are categorical
  predictors handled properly yet.  Our version does handle multiple
  response variables, however.  As it is not well-tested, we would like
  to hear of any bugs.}
\author{
  Trevor Hastie and Robert Tibshirani
}
\references{
  J. Friedman, ``Multivariate Additive Regression Splines''.
  \emph{Annals of Statistics}, 1991.
}
\seealso{
  \code{\link{predict.mars}},
  \code{\link{model.matrix.mars}}
}
\examples{
data(trees)
fit1 <- mars(trees[,-3], trees[3])
showcuts <- function(obj)
{
  tmp <- obj$cuts[obj$sel, ]
  dimnames(tmp) <- list(NULL, names(trees)[-3])
  tmp
}
showcuts(fit1)

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
