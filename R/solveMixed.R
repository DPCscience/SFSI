
# X=Z=K=indexK=h2=NULL; BLUE=TRUE; BLUP=TRUE;
# method="ML";  return.Hinv = FALSE; tol=1E-5; maxIter=1000; interval=c(1E-10,1E10)

solveMixed <- function(y, X = NULL, Z = NULL, K = NULL, indexK = NULL,
                      h2 = NULL, BLUE = TRUE, BLUP = TRUE, method = "ML",
                      return.Hinv = FALSE, tol = 1E-5, maxIter = 1000,
                      interval = c(1E-10,1E10))
{
  method <- match.arg(method)
  if((BLUE+BLUP)==0) stop("Either 'BLUE' or 'BLUP', or both must be 'TRUE'")
  if(!is.vector(y)) stop("Object 'y' must be a vector\n")

  if(is.character(K)){
    K <- readBinary(K,indexRow=indexK,indexCol=indexK)
  }

  indexOK <- which(!is.na(y))
  anyNA <- any(is.na(y))

  if(is.null(X))
  {
    X <- stats::model.matrix(~1,data=data.frame(rep(1,length(y))))
  }else{
    if(is.vector(X)){
      X <- stats::model.matrix(~X)
      if(ncol(X)>2)  colnames(X)[-1] <- substr(colnames(X)[-1],2,nchar(colnames(X)[-1]))
    }
  }
  stopifnot(nrow(X) == length(y))

  if(!is.null(Z))
  {
    if(!is.matrix(Z)) stop("Object 'Z' must be a matrix with nrow(Z)=n and ncol(Z)=nrow(K)\n")
    K <- Z%*%K%*%t(Z)
  }
  stopifnot(nrow(K) == length(y))
  stopifnot(ncol(K) == length(y))

  if(anyNA){
    K12 <- K[indexOK,-indexOK,drop=FALSE]
  }else K12 <- NULL

  out <- eigen(K[indexOK,indexOK])
  d <- out$values
  U <- out$vectors

  uHat <- rep(NA,length(y))
  if(anyNA){
    X <- X[indexOK, ,drop=FALSE]
    y <- y[indexOK]
  }
  n <- length(y)

  stopifnot(nrow(U) == length(y))
  stopifnot(ncol(U) == length(d))

  Uty <- drop(crossprod(U,y))
  UtX <- crossprod(U,X)

  convergence <- NULL
  if(is.null(h2))
  {
    bb <- exp(seq(log(interval[1]),log(interval[2]),length=15))
    for(i in 1:(length(bb)-1))
    {
      tmp <- try(uniroot(f=dloglik,interval=c(bb[1],bb[i+1]),n=n,
                 Uty=Uty,UtX=UtX,d=d,tol=tol,maxiter=maxIter,trace=2),
              silent = TRUE)
      if(class(tmp) == "list"){
        convergence <- tmp$iter <= maxIter
        lambda0 <- tmp$root
        break
      }
    }
    if(is.null(convergence))
      stop("Values at end points not of opposite sign. Try using function 'solveMixed' with a\n",
           "\tlarger interval or an smaller value of 'tol' or larger 'maxIter' parameters")
  }else lambda0 <- h2/(1-h2)

  dbar <- 1/(lambda0*d+1)
  qq1 <- t(Uty*dbar)%*%UtX
  qq2 <- solve(sweep(t(UtX),2L,dbar,FUN="*")%*%UtX)
  ytPy <- drop(sum(dbar*Uty^2)-qq1%*%qq2%*%t(qq1))

  # Compute BLUE
  if(BLUE){
    bHat <- drop(qq2%*%t(qq1))
    Xb <- as.vector(X%*%bHat)
  }else{
    bHat <- NULL
    Xb <- rep(mean(y),n)
  }

  # Compute BLUP: uHat = AZ' V^{-1} (y-X*b)
  if(BLUP){
    Hinv <- sweep(U,2L,lambda0*dbar,FUN="*") %*% t(U)  # V^{-1}
    uHat[indexOK] <- drop(sweep(U,2L,d*lambda0*dbar,FUN="*")%*%t(U)%*%(y-Xb))
    if(anyNA){
      uHat[-indexOK] <- crossprod(K12,Hinv)%*%(y-Xb)
    }
    varE <- ytPy/n
    varU <- lambda0*varE
    h2 <- varU/(varU + varE)

    if(!is.null(convergence)){
      if(!convergence){
        warning("Convergence was not reached in the 'EMMA' algorithm",immediate.=TRUE)
        varE <- varU <- h2 <- uHat <- NULL
      }
    }
  }else{
    varE <- varU <- h2 <- uHat <- Hinv <- NULL
  }

  if(!return.Hinv) Hinv <- NULL

  # Compute the log-likelihood, equation (3) in Zhou and Stephens, 2012.
  loglikelihood <- n/2*log(n/(2*pi))-n/2-0.5*sum(log(lambda0*d+1))-n/2*log(ytPy)

  out <- list(varE=varE, varU=varU, h2=h2, b=bHat, u=uHat,
              Hinv=Hinv, convergence=convergence, method=method,
              loglikelihood=loglikelihood)
  return(out)
}
