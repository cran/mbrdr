\name{mbrdr}
\alias{mbrdr}
\alias{mbrdr.compute}

\title{Main function for model-based response dimension reduction regression}
\description{
 This is the main function in the mbrdr package.  It creates objects of class
 mbrdr to estimate the response mean subspace and perform tests concerning
 its dimension.  Several helper functions that require a mbrdr object can then
 be applied to the output from this function.
}
\usage{
mbrdr (formula, data, subset, na.action = na.fail, weights, ...)

mbrdr.compute (y, x, weights, method = "upfrr", ...)
 }

\arguments{
 \item{formula}{a two-sided formula like \code{cbind(y1,y2,y3,y4)~x1+x2+x3}, where the left-side
 variables are a matrix of the response variables, and the right-hand side
 variables represent the predictors.  The left-hand side of the formula must be a matrix,
since the package reduces the dimension of the responses variables.
}
 \item{data}{ an optional data frame containing the variables in the model.
           By default the variables are taken from the environment from
           which `mbrdr' is called.}
 \item{subset}{an optional vector specifying a subset of observations to be
          used in the fitting process.}
 \item{weights}{an optional vector of weights to be used where appropriate.  In the
          context of dimension reduction methods, weights are used to obtain
          elliptical symmetry, not constant variance.
          %See \code{\link{mbrdr.weights}}.
          }
 \item{na.action}{a function which indicates what should happen when the data
          contain `NA's.  The default is `na.fail,' which will stop calculations.
          The option 'na.omit' is also permitted, but it may not work correctly when
          weights are used.}
 \item{x}{The design matrix.  This will be computed from the formula by \code{dr} and then
 passed to \code{dr.compute}, or you can create it yourself.}
 \item{y}{The response vector or matrix}
 \item{method}{This character string specifies the method of fitting.
The default is \code{"upfrr"}. The options include \code{"yc"}, \code{"prr"}, \code{"pfrr"}.
 Each method may have its own additional arguments, or its  own defaults; see the details below for more information.}
 \item{\dots}{For \code{mbrdr}, all additional arguments passed to \code{mbrdr.compute}.
For  \code{mbrdr.compute}, additional arguments may be required for particular dimension reduction method.
For  example,  \code{numdir} is the maximum number of directions to compute, with default equal to 4.
Other methods may have other defaults.}
}

\details{
The general regression problem mainly focuses on studying \eqn{E(y|x)},
the conditional mean of a response \eqn{y} given a set of predictors \eqn{x},
where y is \eqn{r}-dimensional response variables with \eqn{r geq 2} and

This function provides methods for estimating the response dimension subspace of a general regression problem.
That is, we want to find a \eqn{r \times d}{r by d} matrix \eqn{B} of minimal rank \eqn{d} such that
\deqn{E(y|x)=E(P(B)y|x)}, where P(B) is an orthogonal projections onto the column space of B.
Both the dimension \eqn{d} and the subspace \eqn{P(B)} are unknown.
These methods make few assumptions.

For the methods \code{"yc"}, \code{"prr"}, \code{"pfrr"} and
\code{"upfrr"},  \eqn{B} is estimated and returned.
And, only for \code{"pfrr"} and \code{"upfrr"},
chi-squared test results for estimating \eqn{d} is provided.

Weights can be used, essentially to specify the relative
frequency of each case in the data.

The option \code{fx.choice} is required to fit \code{"pfrr"} and \code{"upfrr"} and has the following four values.

\code{fx.choice=1}: This is default and returns the original predictor matrice X, centered at zero as fx.

\code{fx.choice=2}: This returns the original predictor matrice X, centered at zero and its squared values.

\code{fx.choice=3}: This returns the original predictor matrice X, centered at zero and its exponentiated values.

\code{fx.choice=4}: This clusters X with K-means algoritm with the number of clusters equal to the value in \code{nclust}.
Then, the cluster results are expanded to \eqn{\code{nclust}-1} dummy variables, like factor used in \code{lm} function. Finally, it returns \code{nclust-1} categorical basis. The option of \code{nclust} works only with \code{fx.choice=4}.

}

\value{
mbrdr returns an object that inherits from mbrdr (the name of the type is the
value of the \code{method} argument), with attributes:
  \item{y}{The response matrix}
  \item{x}{The design matrix}
  \item{weights}{The weights used, normalized to add to n.}
  \item{cases}{Number of cases used.}
  \item{call}{The initial call to \code{mbrdr}.}
  \item{evectors}{The eigenvectors from kernel matrices to estimate \eqn{B} computed from each response dimension reduction methods. It is the estimate of \eqn{B}.}
  \item{evalues}{The eigenvalues corresponding to the eigenvectors.}
  \item{stats}{This is the dimension test statistics for \eqn{pfrr} and \code{"upfrr"}. It is the cumulatative sum of the eigenvalues for \code{"yc"} and \code{"prr"}}
  \item{fx}{This returns the user-selection of fx for \code{"pfrr"} and \code{"upfrr"}.}
  \item{numdir}{The maximum number of directions to be found.  The output
value of numdir may be smaller than the input value.}
   \item{method}{the dimension reduction method used.}
}

\references{

Yoo, JK. (2018). Response dimension reduction: model-based approach.
\emph{Statistics : A Journal of Theoretical and Applied Statistic}, 52, 409-425. \code{"prr"} and \code{"pfrr"}

Yoo, JK. (2019). Unstructured principal fitted response reduction in multivariate regression.
\emph{ Journal of the Korean Statistical Society}, 48, 561-567.  \code{"upfrr"}

Yoo, JK. and Cook, R. D.  (2008), Response dimension reduction for the conditional mean in multivariate regression.
\emph{Statistics and Probability Letters}, 47, 381-389.  \code{"yc"}.
}

\author{Jae Keun Yoo, <peter.yoo@ewhat.ac.kr>.} %



\examples{
data(mps)
# default fitting method is "upfrr"
s0 <- mbrdr(cbind(A4, B4, A6, B6)~AFDC+Attend+B+Enrol+HS+Minority+Mobility+Poverty+PTR, data=mps)
summary(s0)

# Refit, using different choice of fx.
summary(s1 <- update(s0, fx.choice=2))

# Refit again, using pfrr with fx.choice=2
summary(s2<-update(s1, method="pfrr", fx.choice=1))

# Refit, using prr, which does not require the choice of fx.
summary(s3<- update(s1,method="prr"))

# fit using Yoo-Cook method:
summary(s4 <- update(s1,method="yc"))
}

\keyword{regression }%-- one or more ...
