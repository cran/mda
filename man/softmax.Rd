\name{softmax}
\alias{softmax}
\title{
identify the maximum in each row of a matrix
}
\usage{
softmax(x, gap=F)
}
\arguments{
\item{x}{
a matrix of mode numeric
}
\item{gap}{
if \code{TRUE}, the difference between the largest and next largest column is returned.
}}
\value{
a factor with levels the column labels of \code{x} and values the columns corresponding to the maximum column. If \code{gap = TRUE} a list is returned, the second component of which is the difference between the largest and next largest column of \code{x}
}
\seealso{
\code{\link{predict.fda}}, \code{\link{confusion}}, \code{\link{fda}} 
}
% Converted by Sd2Rd version 0.2-a5.
