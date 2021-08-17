\name{mbrdr.x}
\alias{mbrdr.x}
\alias{mbrdr.y}

\title{ Accessor functions for data in dr objects }
\description{
  Accessor functions for dr objects.
}
\usage{
mbrdr.x(object)
mbrdr.y(object)
}

\arguments{
  \item{object}{ An object that inherits from \code{mbrdr}. }
}

\value{
Returns a component of a dr object.  \code{mbrdr.x} returns the predictor
matrix reduced to full rank by dropping trailing columns; \code{mbrdr.y}
returns the response vector/matrix.
}

\author{ Jae Keun Yoo, <peter.yoo@ewha.ac.kr> }

\seealso{  \code{\link{mbrdr}}. }

\keyword{internal}% at least one, from doc/KEYWORDS
