\name{polyreg}
\alias{polyreg}
\title{
simple minded polynomial regression
}
\usage{
polyreg(x, y, w, degree=1, monomial=F, \dots)


}
\arguments{
\item{x}{
predictor matrix
}
\item{y}{
response matrix
}
\item{w}{
optional (positive) weights
}
\item{degree}{
total degree of polynomial basis (default is 1)
}
\item{monomial}{
if \code{TRUE} a monomial basis is used (no cross terms). default is \code{FALSE}
}}
\value{
a polynomial regression fit, containing the essential ingredients
for its predict method.
}
\seealso{
\code{\link{predict.polyreg}}
}
% Converted by Sd2Rd version 0.2-a5.
