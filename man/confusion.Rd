\name{confusion}
\alias{confusion}
\title{
compute the confusion matrix between two factors, or for an fda or mda object
}
\usage{
confusion(predict, true)
}
\arguments{
\item{predict}{
the predicted factor, or an fda or mda model
}
\item{true}{
the true factor, or else (in the case where predict is a model object) a data frame (list) containing the test data.
}}
\value{
essential table(predict, true), but with some useful attribute(s)
}
\seealso{
\code{\link{fda}}, \code{\link{predict.fda}}
}
\examples{
data(iris)
irisfit <- fda(Species ~ Sepal.Length + Sepal.Width + Petal.Length,
    data = iris, method=mars)
confusion(predict(irisfit, iris), iris[,"Species"])
}
% Converted by Sd2Rd version 0.2-a5.
