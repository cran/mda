\name{coef.fda}
\alias{coef.fda}
\alias{coef.mda}
\title{
Produce coefficients for an fda or mda object
}
\description{
a method for coef for extracting the canonical coefficients from an fda or mda object
}
\usage{
\method{coef}{fda}(object, \dots)
}
\arguments{
  \item{object}{an \code{fda} or \code{mda} object.
}
  \item{\dots}{
not relevant
}
}
\details{
See the references for details.
}
\value{
A coefficient matrix}
\author{
  Trevor Hastie and Robert Tibshirani
}
\seealso{
  \code{\link{predict.fda}},
  \code{\link{plot.fda}},
  \code{\link{mars}},
  \code{\link{bruto}},
  \code{\link{polyreg}},
  \code{\link{softmax}},
  \code{\link{confusion}},
  %% \code{\link{coefficients.fda}}
}
\references{
  ``Flexible Disriminant Analysis by Optimal Scoring''  by Hastie,
  Tibshirani and Buja, 1994, JASA, 1255-1270.

  ``Penalized Discriminant Analysis'' by Hastie, Buja and Tibshirani, 1995,
  Annals of Statistics, 73-102.

  ``Elements of Statisical Learning - Data Mining, Inference and
  Prediction'' (2nd edition, Chapter 12) by Hastie, Tibshirani and
  Friedman, 2009, Springer
}
\examples{
data(iris)
irisfit <- fda(Species ~ ., data = iris)
coef(irisfit)
mfit=mda(Species~.,data=iris,subclass=2)
coef(mfit)
}\keyword{classif}

