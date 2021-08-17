\name{choose.fx}
\alias{choose.fx}


\title{choose fx for principal fitted response reduction and unstructured principal fitted response reduction }

\description{
 Returns a \eqn{n \times q}  matrix used in  principal fitted response reduction and unstructured principal fitted response reduction.}

\usage{
choose.fx(X, fx.choice=1, nclust = 5)
}

\arguments{
  \item{X}{  \eqn{n \times p} predictor matrix }
  \item{fx.choice}{four choices for fx; see below }
  \item{nclust}{the number of clusters; see below }
}
\details{
Both of principal fitted response reduction and unstructured principal fitted response reduction require
a choice of fx. The function will return one of four choices of fx, which are popular candidates among many.

\code{fx.choice=1}: This is default and returns the original predictor matrice X, centered at zero as fx.

\code{fx.choice=2}: This returns the original predictor matrice X, centered at zero and its squared values.

\code{fx.choice=3}: This returns the original predictor matrice X, centered at zero and its exponentiated values.

\code{fx.choice=4}: This clusters X with K-means algoritm with the number of clusters equal to the value in \code{nclust}.
Then, the cluster results are expanded to \eqn{\code{nclust}-1} dummy variables, like factor used in \code{lm} function.
Finally, it returns \code{nclust-1} categorical basis. The option of \code{nclust} works only with \code{fx.choice=4}.
}
\value{
 A \eqn{n \times q}  matrix for fx.
}
\author{Jae Keun Yoo, peter.yoo@ewha.ac.kr}

\examples{
data(mps)
X <- mps[,c(5:6,8:14)]
choose.fx(X)

choose.fx(X, fx.choice=2)

choose.fx(X, fx.choice=4, nclust=3)
}

\keyword{ internal }
