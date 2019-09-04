#' Least Angle Regression to solve the LASSO-type problem
#'
#' Computes the entire LASSO solution for the regression coefficients, starting from zero, to the
#' least squares estimates, via the Least Angle Regression (LARS) algorithm (Efron, 2004). It uses as inputs
#' a variance matrix among predictors and a covariance vector between response and predictors.
#' 
#' The regression coefficients \eqn{\beta=(\beta_1,...,\beta_p)'} are estimated as function of the variance matrix among
#' predictors (\eqn{XtX}) and the covariance vector between response and predictors (\eqn{Xty}) using 'covariance updates'
#' to solve the optimization function
#' \deqn{-Xty' \beta + 1/2\beta' (XtX)\beta + 1/2\lambda||\beta||_2^2}
#' where \eqn{\lambda} is the penalization parameter
#' @return  List with the following elements:
#' \itemize{
#'   \item \code{beta}: vector of regression coefficients.
#'   \item \code{lambda}: penalty of LASSO-type problem for all the sequence of coefficients.
#'   \item \code{df}: degrees of freedom, number of non-zero predictors at each solution.
#'   \item \code{sdx}: vector of standard deviation of predictors.
#' }
#' @param XtX Variance-covariance matrix among predictors
#' @param Xty Covariance vector between response variable and predictors
#' @param method One of:
#' \itemize{
#'  \item \code{'LAR'}: Computes the entire sequence of all coefficients. Values of lambdas are calculated at each step.
#'  \item \code{'LAR-LASSO'}: Similar to \code{'LAR'} but solutions when a predictor leaves the solution are also returned.
#' }
#' Default is \code{method='LAR'}
#' @param eps An effective zero. Default is the machine precision
#' @param maxDF Maximum number of predictors in the last lars solution.
#' Default \code{maxDF=NULL} will calculate solution for all the predictors
#' @param scale \code{TRUE} or \code{FALSE} to recalculate the matrix \code{XtX} for variables with unit variance 
#' (see \code{help(scale_crossprod)}) and scale \code{Xty} by the standard deviation of the corresponding predictor
#' taken from the diagonal of \code{XtX}
#' @param verbose \code{TRUE} or \code{FALSE} to whether printing each lars step
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate variables
#' n = 500; p=200;  rho=0.65
#' X = matrix(rnorm(n*p),ncol=p)
#' eta = scale(X%*%rnorm(p))  # signal
#' e =  rnorm(n)              # noise
#' y = rho*eta + sqrt(1-rho^2)*e
#'
#' # Training and testing sets
#' pTST = 0.3      # percentage to predict
#' tst = sample(1:n,floor(pTST*n))
#' trn = (1:n)[-tst]
#'
#' # Calculate covariances in training set
#' P = var(X[trn,])
#' rhs = as.vector(cov(y[trn],X[trn,]))
#'
#' # Run the penalized regression
#' fm = lars2(P,rhs,method="LAR-LASSO",verbose=TRUE)
#'
#' # Regression coefficients
#' beta = as.matrix(fm$beta)
#'
#' # Predicted values in training and testing set
#' yHat_TRN =  X[trn,] %*% t(beta)
#' yHat_TST =  X[tst,] %*% t(beta)
#'
#' par(mfrow=c(1,2))
#' plot(fm$df,cor(y[trn],yHat_TRN)[1,],main="Training set")
#' plot(fm$df,cor(y[tst],yHat_TST)[1,],main="Testing set")
#' @export
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos. Adapted from 'lars' package (Hastie & Efron, 2013)
#' @references
#' \itemize{
#' \item \insertRef{Efron2004}{SFSI}
#' \item \insertRef{Friedman2010}{SFSI}
#' \item \insertRef{Hastie2013}{SFSI}
#' \item \insertRef{Tibshirani1996}{SFSI}
#' }
#' @keywords lars2

lars2 <- function(XtX, Xty, method=c("LAR","LAR-LASSO"), maxDF=NULL,
    eps=.Machine$double.eps,scale=TRUE,verbose=FALSE)
{
  method <- match.arg(method)

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
    XtX <- scale_crossprod(XtX)
    sdx <- XtX$sdX
    XtX <- XtX$XtX
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
  return(out)
}
