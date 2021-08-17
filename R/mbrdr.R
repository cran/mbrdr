.packageName <- "mbrdr"
#####################################################################
#     mbrdr is the primary function
#####################################################################

mbrdr <- function(formula, data, subset, na.action=na.fail, weights, ...){
    mf <- match.call(expand.dots=FALSE)
    mf$na.action <- na.action
    mf$... <- NULL # ignore ...
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")

    Y <- model.extract(mf,"response")
    if (dim(Y)[2]<2)     #check multi-dimensionality of Y
        stop("The responses must be multivariate !")

    X <- model.matrix(mt, mf)
    y.name <- if(is.matrix(Y)) colnames(Y) else
        as.character(attr(mt, "variables")[2])
    int <- match("(Intercept)", dimnames(X)[[2]], nomatch=0)
    if (int > 0) X <- X[, -int, drop=FALSE] # drop the intercept from X

    weights <- mf$"(weights)"
    if(is.null(weights)){                   # create weights if not present
         weights <- rep(1,dim(X)[1])} else
        {
         if(any(weights<0))stop("Negative weights")
         pos.weights <- weights > 1.e-30
         weights <- weights[pos.weights]
         weights <- dim(X)[1]*weights/sum(weights)
         X <- X[pos.weights,]
         Y <- if(is.matrix(Y)) Y[pos.weights,] else Y[pos.weights]}

    ans <- mbrdr.compute(Y, X, weights=weights,...)
    ans$call <- match.call()
    ans$y.name <- y.name
    ans
}


mbrdr.compute <- function(y, x, weights, method="upfrr", ...){
    if (NROW(y) != nrow(x))     #check lengths
        stop("The response and predictors have differing number of observations")

    classname<- c(method)
    genclassname<-"mbrdr"
    sweights <- sqrt(weights)

#initialize the object and then call the fit method
    ans <- list(y=y, x=x, weights=weights,method=method,cases=NROW(y))
    class(ans) <-  c(classname,genclassname)   #set the class
    ans <- mbrdr.fit(object=ans, ...)   # ... contains args for the method
    class(ans) <-  c(classname,genclassname) # reassign class name
    ans #return the object
}


#######################################################
#    accessor functions
#######################################################

mbrdr.x <- function(object) {object$x}
mbrdr.wts <- function(object) {object$weights}
mbrdr.y.name <- function(object) {object$y.name}

#####################################################################
#     Fitting function
#####################################################################
mbrdr.fit <- function(object, numdir=4,  fx.choice=1, fx=NULL, nclust=5,...) UseMethod("mbrdr.fit")

mbrdr.fit.default <-function(object, numdir=4, fx.choice=1,fx=NULL, nclust=5,...){
    M <-mbrdr.M(object, ...)  # Get the kernel matrix M, method specific
    evalues <- M$evalues
    evectors <- as.matrix(M$M)
    numdir <- min( M$numdir, dim(evectors)[2] )
    stats <- M$stats
    if (is.null(M$fx)) {fx <- NULL}
      else{ fx <- M$fx}
    names <- colnames(mbrdr.y(object))[1:dim(evectors)[1]]
    dimnames(evectors)<-
         list(names, paste("Dir", 1:NCOL(evectors), sep=""))
    aa<-c( object, list(evectors=evectors, evalues=evalues, stats=stats, fx=fx, numdir=numdir))
    class(aa) <- class(object)
    return(aa)
}


#####################################################################
###
###  mbrdr methods. Each method REQUIRES a dr.M function, and
###  may also have a dr.y function and an function method.
###
#####################################################################

mbrdr.M <- function(object,   ...){UseMethod("mbrdr.M")}
mbrdr.y <- function(object) {UseMethod("mbrdr.y")}
mbrdr.y.default <- function(object){object$y}
mbrdr.test <- function(object, nd){  UseMethod("mbrdr.test")}
mbrdr.test.default <-function(object, nd) {NULL}


#####################################################################
#     Yoo and Cook (2008) method
#####################################################################
mbrdr.M.yc <- function(object, num.d=NULL,  ...) {

 y <- mbrdr.y(object)
 x <- mbrdr.x(object)

 n <- dim(y)[1]
 r <- dim(y)[2]

 if (is.null(num.d) ) num.d <- (r-1)

 sigmay=cov(y)
 sigmax=cov(x)

 Sigmahat = solve(sigmay) %*% cov(y,x) %*% solve(sigmax)

 eg_S<- eigen(Sigmahat%*%t(Sigmahat));   ev_S <- eg_S$vectors
 evectors <- as.matrix(ev_S[,1:num.d])
 evalues <- eg_S$values[1:num.d]
 cumsum.evalues <- rev(cumsum(rev(evalues)))
 ans <- list(M=evectors, stats=n*cumsum.evalues, evalues=evalues, numdir=num.d)
 return(ans)
}


#####################################################################
#     Principal Response Reduction (PRR)
#####################################################################
mbrdr.M.prr <- function(object, num.d=NULL,...) {

 y <- mbrdr.y(object)
 n <- dim(y)[1]
 r <- dim(y)[2]

 if (is.null(num.d) ) num.d <- (r-1)

 Sigmahat = cov(y)
 eg_S<- eigen(Sigmahat);   ev_S <- eg_S$vectors
 evectors <- as.matrix(ev_S[,1:num.d])
 evalues <- eg_S$values[1:num.d]
 cumsum.evalues <- rev(cumsum(rev(evalues)))
 ans <- list(M=evectors, stats=n*cumsum.evalues, evalues=evalues, numdir=num.d)
 return(ans)
}


#####################################################################
#     Principal Fitted Response Reduction (PFRR)
#####################################################################
mbrdr.M.pfrr <- function(object, num.d=NULL, fx.choice=1, nclust=5, fx=NULL,...) {

####################
## Inner functions begin
####################
  lik <- function(h, S, S_res, n) {
    p <- dim(as.matrix(h))[1]; d <- dim(as.matrix(h))[2]
    h <- eigen(h %*% t(h))$vectors[, 1:d]
    h0 <- eigen(h %*% t(h))$vectors[, (d+1):p]
    lik <- (-n/2)*log( det( t(h0) %*% S %*% h0 ) ) + (-n/2)*log( det( t(h) %*% S_res %*% h ) )
    return(lik)}

  seq.dir <- function(cand_mat, d=1, S, S_res, n){
   h <- h1 <- NULL
   r <- dim(cand_mat)[2]
   f.lik <- -(n/2)*log(det(S_res))
   for (i in 1:d){
     l.lik <- NULL
      for(j in 1:r){
        l.lik[j] <- lik(cbind(h1, cand_mat[,j]), S, S_res, n)
      }
    stat <- 2*(f.lik - l.lik)
    sel <- which(stat >= 0 ); l.lik <- l.lik[sel]; stat <- stat[sel]
    cand_mat <- cand_mat[, sel]
    h <- cbind(h,cand_mat[,which.max(l.lik)])
    h1 <- h
    cand_mat <- cand_mat[, -(which.max(l.lik))]
    r <- length(sel)-1
    }
    out <- list(dir=h, loglik=max(l.lik), stat=min(stat) )
    return(out)}
####################
## Inner functions end
####################

 y <- mbrdr.y(object)
 x <- mbrdr.x(object)

 n <- dim(y)[1]
 r <- dim(y)[2]

 if (is.null(num.d) ) num.d <- (r-1)

 if (fx.choice>4 & is.null(fx)) stop("fx.choice must be less than or equal to 4. If you want to other candidates for fx, use the fx option") else {
      if (!is.null(fx)) fx=fx else fx <- choose.fx(x, fx.choice=fx.choice, nclust=nclust)}

 Sigmas<-SIGMAS(y, fx)
 Sigmahat <- Sigmas$Sigmahat
 Sigmahat_fit<-Sigmas$Sigmahat_fit
 Sigmahat_res<-Sigmas$Sigmahat_res

 eg_S<- eigen(Sigmahat);   ev_S <- eg_S$vectors
 eg_fit <- eigen(Sigmahat_fit);  ev_fit <- eg_fit$vectors
 eg_res <- eigen(Sigmahat_res);  ev_res <- eg_res$vectors

 cand <- cbind(ev_fit, ev_S, ev_res)

 logliks <- NULL
 for (i in 1:num.d ){
   logliks[i] <- seq.dir(cand, d=i, S=Sigmahat, S_res=Sigmahat_res, n=n)$loglik
  }

 evectors <- seq.dir(cand, d=num.d, S=Sigmahat, S_res=Sigmahat_res, n=n)$dir
 evectors <- qr.Q( qr(evectors) )

 evalues <- c((-(n/2)*log(det(Sigmahat))), logliks)

 full.loglik <- (-(n/2)*log(det(Sigmahat_res)))
 stat <- 2*(full.loglik  - evalues )



 ans <- list(M=evectors, stats=stat, evalues=evalues, fx=fx, numdir=num.d)
 return(ans)
}

###################################################################
#  LRT functions
###################################################################

mbrdr.test.pfrr <-function(object, nd) {
    stat <- object$stats
    y <- mbrdr.y(object)
    fx <- object$fx
    r <- dim(y)[2];     q<- dim(fx)[2]

    st<-df<-pv<-0
    nt <- nd
    for (i in 0:nt-1)
      {st[i+1]<-stat[i+1]
       df[i+1]<-(r-i)*q
       pv[i+1]<-1-pchisq(st[i+1],df[i+1])
      }
    z<-data.frame(cbind(st,df, round(pv,4)))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    dimnames(z)<-list(rr,c("Stat","df","p-value"))
    z
}


#####################################################################
#     UnStructured Principal Fitted Response Reduction (UPFRR)
#####################################################################

mbrdr.M.upfrr <- function(object, num.d=NULL, fx.choice=1, nclust=5, fx=NULL, epsilon=1e-7,...) {

#####################################################################
# UPFRR function
#####################################################################
 upfrr = function(Y, d=d, fx=fx, epsilon=1e-7,...){

 ##** Get dimensions
   n <- dim(Y)[1]
   r <- dim(Y)[2]
   q <- dim(fx)[2]

 ##** Compute the sample mean
   muhat <- array(colMeans(Y),dim=c(1,r))

 ##** Compute centered data
   Y.centered <- as.matrix(Y-array(rep(1,times=n),dim=c(n,1))%*%muhat );

 ##** Sample covariance matrices and derived matrices;
   Sigmas<-SIGMAS(Y, fx); Sigmahat <- Sigmas$Sigmahat;
   Sigmahat_fit<-Sigmas$Sigmahat_fit;Sigmahat_res<-Sigmas$Sigmahat_res;

  ##** Unstructured PFRR

   sqrt_Sigmahat_res<-matpower(Sigmahat_res,0.5); Inv_Sqrt_Sigmahat_res<-solve(sqrt_Sigmahat_res);
   lf_matrix <- Inv_Sqrt_Sigmahat_res%*% Sigmahat_fit %*% Inv_Sqrt_Sigmahat_res;
   Vd<-eigen(lf_matrix, symmetric=TRUE)$vectors[, 1:d];
   Gammahat<-(sqrt_Sigmahat_res)%*%Vd%*%solve( matpower( (t(Vd)%*%Sigmahat_res%*%Vd), 0.5)  ) ;

   ##** Estimate Parameters **#;
    Khat<-diag(0, r);
    if ( d<min(r,q) ) {diag(Khat)[(d+1):min(r,q)]<- eigen(lf_matrix, symmetric=TRUE)$values[(d+1):min(r,q)] };
    if ((r<d) | (r <d)) stop("d needs to be less than min(r, q)");
    Vhat<-eigen(lf_matrix, symmetric=TRUE)$vectors;
    Deltahat<-Sigmahat_res + sqrt_Sigmahat_res  %*% Vhat %*% Khat %*% t(Vhat) %*% sqrt_Sigmahat_res;
    Betahat <- solve(t(Gammahat)%*%solve(Deltahat)%*%Gammahat)%*%t(Gammahat) %*% solve(Deltahat)%*% t(Y.centered) %*% fx%*% solve(t(fx)%*%fx)
    Rhat <- t(Gammahat)%*% solve(Deltahat)%*%t(Y.centered); ### Sufficient statistic

  ### Computing Log-Likelihood
  # fully maximized log-likelihood
    temp0<- -(n*r/2)*(1 + log(2*pi) );
    temp1<- -(n/2)*sum( log( eigen(Sigmahat_res, symmetric=TRUE)$values) );
    temp2<-0;
    if ( d<min(r,q) ) {
	  temp2<- -(n/2)*sum(  log( 1+ eigen(lf_matrix, symmetric=TRUE)$values[(d+1):min(r,q)] )  );
  	}
    loglik=temp0 + temp1 + temp2;
  return(list(Gammahat=Gammahat, loglik=loglik, muhat=muhat, Betahat=Betahat, Deltahat=Deltahat, Rhat=Rhat))
  }
#####################################################################

 y <- mbrdr.y(object)
 x <- mbrdr.x(object)

 n <- dim(y)[1]
 r <- dim(y)[2]

 if (is.null(num.d) ) num.d <- (r-1)

 if (fx.choice>4 & is.null(fx)) stop("fx.choice must be less than or equal to 4. If you want to other candidates for fx, use the fx option") else {
      if (!is.null(fx)) fx=fx else fx <- choose.fx(x, fx.choice=fx.choice, nclust=nclust)}

 logliks <- NULL
 Gammahats <- list(NULL)
 for (i in 1:num.d) {
   temp.fit <- upfrr(y, d=i, fx=fx)
   logliks[i] <- temp.fit$loglik
   Gammahats[[i]] <- temp.fit$Gammahat
 }
  d0.loglik <- upfrr(y, d=0, fx=fx)$loglik
  full.loglik <- upfrr(y, d=r, fx=fx)$loglik
  logliks <- c(d0.loglik, logliks)
  stat <- 2*(full.loglik - logliks)

  ans <- list(M=Gammahats[[num.d]], stats=stat, evalues=logliks, Gammahats=Gammahats, fx=fx, numdir=num.d)
return(ans)
}

###################################################################
#  LRT functions
###################################################################

mbrdr.test.upfrr <-function(object, nd) {
    stat <- object$stats
    y <- mbrdr.y(object)
    fx <- object$fx
    r <- dim(y)[2];     q<- dim(fx)[2]

    st<-df<-pv<-0
    nt <- nd

    for (i in 0:nt-1)
      {st[i+1]<-stat[i+1]
       df[i+1]<-(r-i)*(q-i)
       pv[i+1]<-1-pchisq(st[i+1],df[i+1])
      }
    z<-data.frame(cbind(st,df, round(pv,4)))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    dimnames(z)<-list(rr,c("Stat","df","p-value"))
    z
}


###################################################################
#
#  basic print method for dimension reduction
#
###################################################################
"print.mbrdr" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Eigenvectors:\n")
    evectors<-as.matrix(x$evectors)
    print.default(evectors)
    stats <- t(as.matrix(round(x$stats,4)))
    nt <- length(x$stats)
    rr <-  paste("H0: d=", 0:(nt-1),sep="")
    colnames(stats)<-rr
    rownames(stats) <- c("")
    cat("Statistics:\n")
    print(stats)
    cat("\n")
    invisible(x)
}

###################################################################
#  basic summary method for dimension reduction
###################################################################
"summary.mbrdr" <- function (object, ...){
    z <- object
    ans <- z[c("call")]

    nd <- z$numdir
    stats <- t(as.matrix(round(z$stats,4)))
    nt <- length(z$stats)
    rr <-  paste("H0: d=", 0:(nt-1),sep="")
    colnames(stats)<-rr
    rownames(stats) <- c("")

    ans$evectors <- z$evectors[,1:nd]
    ans$method <- z$method
    ans$weights <- mbrdr.wts(z)
    sw <- sqrt(ans$weights)
    y <- z$y
    ans$n <- z$cases #NROW(z$model)

    ans$test <- mbrdr.test(z, nt)
    ans$stats <- stats
    class(ans) <- "summary.mbrdr"
    ans
}

###################################################################
#
# basic print.summary method for dimension reduction
#
###################################################################
"print.summary.mbrdr" <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")#S: ' ' instead of '\n'
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")

    cat("Method:\n")#S: ' ' instead of '\n'
       cat(paste(x$method, ", n = ", x$n,sep=""))
       if(diff(range(x$weights)) > 0)
                  cat(", using weights.\n") else cat(".\n")
    cat("\n")
    cat("\nEigenvectors:\n")
    print(x$evectors, digits=digits)
    cat("\n")
    cat("\nStatistics:\n")
    print(x$stats,digits=digits)
    cat("\n")
    if (!is.null(x$test)){
      cat("\nAsymp. Chi-square tests for dimension:\n")
      print(as.matrix(x$test),digits=digits)}
    invisible(x)
}

###################################################################
# Required functions
###################################################################

## function 1
choose.fx = function(X, fx.choice=1, nclust=5){
  CategoricalBasis<-function(y){
   nobs=length(y); bin.y<-unique(sort(y));
   r<- length(unique(sort(y)))-1; fy<-array(rep(0), c(r, nobs));
   ny<-array(rep(0), c(r, 1)); # Use to count of points falling into the bins, n.y[i] is for bin i
    for (i in 1:r){ fy[i,]<-sapply(y, function(x) (x==bin.y[i]) ) }
   return(fy);
  }
  Xc <- scale(X, scale=F)
  if (fx.choice==1) {fx = Xc}
  if (fx.choice==2) {fx = cbind(Xc, Xc^2)}
  if (fx.choice==3) {fx = cbind(Xc, exp(Xc))}
  if (fx.choice==4) {slice = kmeans(X, nclust)$clust; fx=t(CategoricalBasis(slice))}
  if (fx.choice >4) {fx= NULL}
return(fx)
}

## function 2
SIGMAS = function(Y, fx)
{

### Get dimensions
if (!is.vector(Y)) {

	nobs <- dim(Y)[1]; npred <- dim(Y)[2];

	mu <- array(colMeans(Y),dim=c(1,npred));
}

if (is.vector(Y)){
	nobs <- length(Y); npred <- 1 ; Y<-array(Y, dim=c(nobs,1));

	mu <- array(mean(Y),dim=c(1,npred));
}

r <- dim(fx)[2]; # length of fx

# Compute centered data
X.centered <- as.matrix( Y-array(rep(1,times=nobs),dim=c(nobs,1))%*%mu ) ;

# Sample covariance matrix
Sigmahat <- 1/nobs*t(X.centered)%*%X.centered;

# Compute the sample covariance matrix of the fitted values
P_F <- fx%*%solve(t(fx)%*%fx)%*%t(fx);

Sigmahat_fit <- 1/nobs*t(X.centered)%*%P_F%*%X.centered; rm(P_F)

Sigmahat_res<-Sigmahat-Sigmahat_fit;

return(list(Sigmahat_fit=Sigmahat_fit, Sigmahat=Sigmahat, Sigmahat_res=Sigmahat_res));

}

## function 3
matpower <-function(M,pow)
{
	if (!is.matrix(M)) stop("The argument is not a matrix")

	if ((dim(M)[1]==1) & (dim(M)[2]==1)) return( as.matrix(M^pow) );

	svdM<-svd(M); return(svdM$u %*% diag(c(svdM$d)^pow) %*% t(svdM$v));
}


