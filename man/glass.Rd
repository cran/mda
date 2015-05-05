\name{glass}
\alias{glass}
\title{Glass Identification Database}
\description{
  The \code{glass} data frame has 214 observations and 10 variables,
  representing glass fragments.
}
\usage{data(glass)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{RI}{refractive index}
    \item{Na}{weight percent in corresponding oxide}
    \item{Mg}{weight percent in corresponding oxide}
    \item{Al}{weight percent in corresponding oxide}
    \item{Si}{weight percent in corresponding oxide}
    \item{K}{weight percent in corresponding oxide}
    \item{Ca}{weight percent in corresponding oxide}
    \item{Ba}{weight percent in corresponding oxide}
    \item{Fe}{weight percent in corresponding oxide}
    \item{Type}{Type of glass:
      \describe{
	\item{1}{building\_windows\_float\_processed,}
	\item{2}{building\_windows\_non\_float\_processed,}
	\item{3}{vehicle\_windows\_float\_processed,}
	\item{4}{vehicle\_windows\_non\_float\_processed (none in this
	  database),}
	\item{5}{containers,}
	\item{6}{tableware,}
	\item{7}{headlamps}
      }
    }
  }
}
\source{
  P. M. Murphy and D. W. Aha (1999),
  UCI Repository of Machine Learning Databases,
  \url{http://ftp.ics.uci.edu/pub/machine-learning-databases}
}
\keyword{datasets}
