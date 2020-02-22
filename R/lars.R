
lars2 <- function(XtX, Xty, method=c("LAR","LAR-LASSO"), maxDF=NULL,
    eps=.Machine$double.eps,scale=TRUE,verbose=FALSE)
{
  method <- match.arg(method)

  Xty <- as.vector(Xty)
  p <- length(Xty)
  if(length(XtX) != p^2)
    stop("Incompatible dimensions between 'XtX' and 'Xty'")
  if((sum(dim(XtX))/2)^2 != length(XtX))
      stop("Object 'XtX' must be a squared matrix. LARS can not be implemented")

  if(!is.double(Xty)) Xty <- as.double(Xty)
  if(!is.double(XtX)) XtX <- apply(XtX,2,as.double)

  im <- inactive <- seq(p)
  Sp <- nchar(p)
  textPrint <- c(" Step","\tSec/Step","\tVariable")

  if(scale){
    sdx <- sqrt(diag(XtX))
    XtX <- scale_cov(XtX)  # Equal to XtX=cov2cor(XtX) but faster
    Xty <- Xty/sdx
  }else{
    sdx <- rep(1,p)
  }

  ignores <- NULL
  if(is.null(maxDF))  maxDF <- p
  beta <- matrix(0,maxDF*8,p)
  lambda <- double(maxDF*8)
  active <- NULL
  drops <- FALSE
  Sign <- NULL
  R <- NULL
  k <- 0
  time <- proc.time()[3]
  while((length(active) < maxDF) & (length(active) < (p-length(ignores))))
  {
    Cov <- Xty[inactive]
    Cmax <- max(abs(Cov))
    if(Cmax < eps*100){
      if(verbose) cat(" Max |corr| = 0; exiting...\n")
      break
    }
    k <- k+1
    lambda[k] <- Cmax
    if(!any(drops))
    {
      new <- abs(Cov) >= Cmax-eps
      Cov <- Cov[!new]
      new <- inactive[new]
      for(inew in new)
      {
        R <- upDateR(XtX[inew,inew],R,drop(XtX[inew,active]),eps=eps)
        if(attr(R,"rank")==length(active))
        {
          nR <- seq(length(active))
          R <- R[nR,nR,drop=FALSE]
          attr(R,"rank") <- length(active)
          ignores <- c(ignores,inew)
          if(verbose){
            cat(" LARS Step ",k,":\t Variable", inew,"\tcollinear; dropped for good\n",sep="")
          }
        }else{
          active <- c(active,inew)
          Sign <- c(Sign,sign(Xty[inew]))
          if(verbose){
            cat("--------------------------------------------------------------\n")
            tmp <- proc.time()[3]
            cat(paste(textPrint,"=",c(sprintf("%*d",Sp,k),sprintf('%.*f',4,tmp-time),
                sprintf("%*d",Sp,inew)))," added\n",sep="")
            time <- tmp
          }
        }
      }
    }
    Gi1 <- backsolve(R,backsolvet(R,Sign))
    A <- 1/sqrt(sum(Gi1*Sign))
    w <- A*Gi1

    if((length(active) >= (p-length(ignores)))){
        gamhat <- Cmax/A
    }else{
      a <- drop(w %*% XtX[active, -c(active,ignores),drop=FALSE])
      gam <- c((Cmax-Cov)/(A-a),(Cmax+Cov)/(A+a))
      gamhat <- min(gam[gam > eps],Cmax/A)
    }
    if(method == "LAR-LASSO")
    {
      dropid <- NULL
      b1 <- beta[k,active]
      z1 <- -b1/w
      zmin <- min(z1[z1 > eps],gamhat)
      if(zmin < gamhat){
          gamhat <- zmin
          drops <- z1 == zmin
      }else drops <- FALSE
    }
    beta[k+1,] <- beta[k,]
    beta[k+1,active] <- beta[k+1,active] + gamhat*w
    Xty <- Xty - gamhat*XtX[,active,drop=FALSE]%*%w
    if(method == "LAR-LASSO" && any(drops))
    {
      dropid <- seq(drops)[drops]
      for(id in rev(dropid))
      {
        if(verbose){
          cat("--------------------------------------------------------------\n")
          tmp <- proc.time()[3]
          cat(paste(textPrint,"=",c(sprintf("%*d",Sp,k+1),sprintf('%.*f',4,tmp-time),
            sprintf("%*d",Sp,active[id])))," dropped\n",sep="")
          time <- tmp
        }
        R <- downDateR(R,id)
      }
      dropid <- active[drops]
      beta[k+1,dropid] <- 0
      active <- active[!drops]
      Sign <- Sign[!drops]
    }
    inactive <- im[-c(active, ignores)]
  }
  beta <- beta[seq(k+1), ,drop = FALSE]
  lambda  <-  c(lambda[seq(k)],0)
  if(scale) beta <- scale(beta,FALSE,sdx)
  df <- apply(beta,1,function(x)sum(abs(x)>0))

  out <- list(method=method,beta=Matrix::Matrix(beta, sparse=TRUE),lambda=lambda,df=df,sdx=sdx)
  class(out) <- "SSI"
  return(out)
}
