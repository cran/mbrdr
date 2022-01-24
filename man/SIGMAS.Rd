\name{SIGMAS}
\alias{SIGMAS}
\title{compute all required SIGMA matrices for \code{"pfrr"} and \code{"upfrr"}}
\description{
 Returns Sigmahat, Sigmahat_fit and Sigmahat_res  for principal fitted response reduction and unstructured principal fitted response reduction using the choice of fx.}
\usage{
SIGMAS(Y, fx)
}

\arguments{
  \item{Y}{  \eqn{n \times r} response matrix }
  \item{fx}{the chosen fx}
}
\details{
Both of principal fitted response reduction and unstructured principal fitted response reduction require
to compute many SIGMAs. The SIGMAs  are as follows:
Sigmahat =  (Y^{T} Y)/n;
Sigmahat_fit = (Y^{T}  P_fx  Y)/n;
Sigmahat_res = Sigmahat - Sigmahat_fit.
}

\value{
A list of Sigmahat, Sigmahat_fit and  Sigmahat_res.
}

\author{Jae Keun Yoo, peter.yoo@ewha.ac.kr}

\examples{
data(mps)
X <- mps[,c(5:6,8:14)]
Y <- mps[,c(1:4)]
fx1 <- choose.fx(X)
fx2 <- choose.fx(X, fx.choice=4, nclust=3)

SIGMAS(Y, fx1)
SIGMAS(Y, fx2)
}

