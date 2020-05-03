
# X = Z = indexK = h2=U=d=NULL; BLUP=TRUE; K=G0; method="ML"
# return.Hinv = FALSE; tol=1E-5; maxIter=1000; interval=c(1E-9,1E9)
solveMixed <- function(y, X = NULL, Z = NULL, K = NULL, U = NULL, d = NULL,
                        indexK = NULL, h2 = NULL, BLUP = TRUE, method = "ML",
                        return.Hinv = FALSE, tol = 1E-5, maxIter = 1000,
                        interval = c(1E-9,1E9))
{
  method <- match.arg(method)
  eps <- .Machine$double.eps

  if(is.character(K)){
    K <- readBinary(K,indexRow=indexK,indexCol=indexK)
  }
  y <- as.vector(y)
  indexNA <- which(is.na(y))
  indexOK <- which(!is.na(y))
  n <- length(indexOK)

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

  if(is.null(U) & is.null(d))
  {
    if(is.null(Z))
    {
      Z <- diag(length(y))
      if(is.null(K)){
          K <- diag(length(y))
      }
      G <- K
    }else{
      if(!is.matrix(Z)) stop("Object 'Z' must be a matrix")
      if(is.null(K)){
        K <- diag(ncol(Z))
        G <- tcrossprod(Z)  # K = ZIZ'
      }else{
        G <- Z%*%K%*%t(Z)
      }
    }
    stopifnot(nrow(G) == length(y))
    stopifnot(ncol(G) == length(y))

    out <- eigen(G[indexOK,indexOK])
    d <- out$values
    U <- out$vectors

  }else{
    if(is.null(U)) stop("You are providing the eigenvalues, but not the eigenvectors")
    if(is.null(d)) stop("You are providing the eigenvectors, but not the eigenvalues")
    if(length(indexNA)>0) stop("No 'NA' values are allowed when passing parameters 'U' and 'd'")
  }

  stopifnot(nrow(U) == n)
  stopifnot(ncol(U) == length(d))

  Uty <- drop(crossprod(U,y[indexOK]))
  UtX <- crossprod(U,X[indexOK, ,drop=FALSE])

  convergence <- lambda0 <- NULL
  if(is.null(h2))
  {
    tmp <- try(uniroot(f=dloglik,interval=interval,n=n,Uty=Uty,UtX=UtX,
                       d=d,tol=tol,maxiter=maxIter,trace=2),
               silent = TRUE)
    if(class(tmp) == "list")
    {
      if(abs(tmp$f.root) <= tol*100){
        convergence <- tmp$iter <= maxIter
        lambda0 <- tmp$root
      }
    }
    if(is.null(convergence))
    {
      # Expand the seeking interval
      bb <- exp(seq(log(eps/10),log(interval[2]^1.7),length=100))
      flag <- TRUE; i <- 1
      while(flag)
      {
        i <- i + 1
        tmp <- try(uniroot(f=dloglik,interval=c(bb[i-1],bb[i]),n=n,Uty=Uty,
                           UtX=UtX,d=d,tol=tol,maxiter=maxIter,trace=4),
                   silent = TRUE)
        if(class(tmp) == "list")
        {
          # If the root is near to zero
          if(abs(tmp$f.root) <= tol*100)
          {
            convergence <- tmp$iter <= maxIter
            if(tmp$root <= interval[1]){
              lambda0 <- interval[1]
            }else{
              if(tmp$root >= interval[2]){
                lambda0 <- interval[2]
              }else lambda0 <-  tmp$root
            }
          }
        }
        #aa <- rep(NA,3)
        #if(class(tmp) == "list") aa=c(tmp$root,tmp$f.root,tmp$estim.prec)
        #cat("Interval ",i,"[",bb[i-1],",",bb[i],"]: root=",aa[1]," f.root=",aa[2]," prec=",aa[3],"\n")
        if(i == length(bb) | !is.null(convergence)) flag <- FALSE
      }

      if(is.null(convergence))
      {  stop("Solution at end points not of opposite sign. Try using function 'solveMixed' with a\n",
               "\tlarger interval or an smaller value of 'tol' or larger 'maxIter' parameters")
      }
    }
  }else lambda0 <- h2/(1-h2)

  dbar <- 1/(lambda0*d+1)
  qq1 <- t(Uty*dbar)%*%UtX
  qq2 <- solve(sweep(t(UtX),2L,dbar,FUN="*")%*%UtX)
  ytPy <- drop(sum(dbar*Uty^2)-qq1%*%qq2%*%t(qq1))

  # Compute BLUE
  bHat <- drop(qq2%*%t(qq1))

  # Compute BLUP: uHat = KZ' V^{-1} (y-X*b)
  uHat <- yHat <- Hinv <- NULL
  if(BLUP)
  {
    yStar <- y[indexOK] - X[indexOK ,,drop=FALSE]%*%bHat
    Hinv <- tcrossprod(sweep(U,2L,lambda0*dbar,FUN="*"),U) # Vinv
    #H <- tcrossprod(sweep(U,2L,d*lambda0*dbar,FUN="*"),U)
    KZt <- tcrossprod(K,Z[indexOK,,drop=FALSE])

    uHat <- drop(KZt%*%Hinv%*%yStar)   # u = K Z' Vinv y*
    yHat[indexOK] <- drop(X[indexOK, ,drop=FALSE]%*%bHat + Z[indexOK, ,drop=FALSE]%*%uHat)
    if(length(indexNA)>0){
      yHat[indexNA] <- drop(X[indexNA,,drop=FALSE]%*%bHat + Z[indexNA, ,drop=FALSE]%*%uHat)
    }

    if(!return.Hinv) Hinv <- NULL

    if(!is.null(convergence))
    { if(!convergence){
        warning("Convergence was not reached in the 'EMMA' algorithm",immediate.=TRUE)
      }
    }
  }

  varE <- ytPy/n
  varU <- lambda0*varE
  h2 <- varU/(varU + varE)

  out <- list(varE = varE, varU = varU, h2 = h2, yHat=yHat, b = bHat, u = uHat,
              Hinv = Hinv, convergence = convergence, method = method)
  return(out)
}
