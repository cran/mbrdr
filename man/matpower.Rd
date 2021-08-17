\name{matpower}
\alias{matpower}
\title{compute the M^power where M is a symmetric matrix.}

\description{
 Returns M^power. }

\usage{
matpower(M, pow)
}

\arguments{
  \item{M}{ symmetric matrix}
  \item{pow}{power}
}
\details{
The function computes M^{power} for a symmetric matrix M. }
\value{
Returns }
\author{Jae Keun Yoo, peter.yoo@ewha.ac.kr}

\examples{
X <- matrix(rnorm(100), c(20,5))
matpower(cov(X), -0.5) ## returns cov(X)^-0.5 %*% cov(X)^-0.5 = cov(X)^-1.
 }

\keyword{ internal }
