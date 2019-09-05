#' Coordinate Descent algorithm to solve the Elastic-Net-type problem
#'
#' Computes the entire Elastic-Net solution for the regression coefficients simultaneously for all
#' values of the penalization parameter using as inputs a variance matrix among predictors and a covariance
#' vector between response and predictors, via the Coordinate Descent (CD) algorithm (Friedman, 2007).
#' 
#' The regression coefficients \eqn{\beta=(\beta_1,...,\beta_p)'} are estimated as function of the variance matrix among
#' predictors (\eqn{XtX}) and the covariance vector between response and predictors (\eqn{Xty}) by optimizing the function
#' \deqn{-Xty' \beta + 1/2\beta' (XtX)\beta + \lambda J(\beta)}
#' where \eqn{\lambda} is the penalization parameter and \eqn{J(\beta)} is a penalty function given by
#' \deqn{1/2(1-\alpha)||\beta||_2^2 + \alpha||\beta||_1}
#' for \eqn{\alpha} between 0 and 1
#' @return  List object containing the elements:
#' \itemize{
#'   \item \code{beta}: vector of regression coefficients.
#'   \item \code{lambda}: sequence of values of lambda used
#'   \item \code{df}: degrees of freedom, number of non-zero predictors at each solution.
#'   \item \code{sdx}: vector of standard deviation of predictors.
#' }
#' @param XtX Variance-covariance matrix among predictors
#' @param Xty Covariance vector between response variable and predictors
#' @param lambda Penalization parameter sequence vector. Default is \code{lambda=NULL}, in this case a decreasing grid of 
#' \code{n='nLambda'} lambdas will be generated starting from a maximum equal to 
#' \tabular{c}{\code{max(abs(Xty)/alpha)}}
#' to a minimum equal to zero. If \code{alpha=0} the grid is generated starting from a maximum equal to 5
#' @param nLambda Number of lambdas generated when \code{lambda=NULL}
#' @param alpha Numeric between 0 and 1 indicating the weights for LASSO (\eqn{\alpha=1}) and Ridge-Regression (\eqn{\alpha=0})
#' @param scale \code{TRUE} or \code{FALSE} to recalculate the matrix \code{XtX} for variables with unit variance 
#' (see \code{help(scale_crossprod)}) and scale \code{Xty} by the standard deviation of the corresponding predictor
#' taken from the diagonal of \code{XtX}
#' @param tol Maximum error between two consecutive solutions of the iterative algorithm to declare convergence
#' @param maxIter Maximum number of iterations to run at each lambda step before convergence is reached
#' @param verbose \code{TRUE} or \code{FALSE} to whether printing each CD step
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate variables
#' n=500; p=200;  rho=0.65
#' X = matrix(rnorm(n*p),ncol=p)
#' signal = rho*scale(X%*%rnorm(p))
#' noise =  sqrt(1-rho^2)*rnorm(n)
#' y = signal + noise
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
#' fm = solveEN(P,rhs,verbose=TRUE)
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
#' @references
#' \itemize{
#' \item \insertRef{Friedman2007}{SFSI}
#' \item \insertRef{Hoerl1970}{SFSI}
#' \item \insertRef{Tibshirani1996}{SFSI}
#' \item \insertRef{Zou2005}{SFSI}
#' }
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#' @keywords solveEN

solveEN <- function(XtX,Xty,lambda=NULL,nLambda=100,alpha=1,scale=TRUE,
  tol=1E-5,maxIter=1000,verbose=FALSE)
{
    p <- length(Xty)
    if(length(XtX) != p^2)
      stop("Incompatible dimensions between 'XtX' and 'Xty'")
    if((sum(dim(XtX))/2)^2 != length(XtX))
      stop("Object 'XtX' must be a squared matrix. The algorithm can not be implemented")

    if(alpha<0 | alpha>1) stop("Parameter 'alpha' must be a number between 0 and 1")

    if(!is.double(Xty)) Xty <- as.double(Xty)
    if(!is.double(XtX)) XtX <- apply(XtX,2,as.double)

    if(scale)
    {
      XtX <- scale_crossprod(XtX)
      sdx <- XtX$sdX
      XtX <- XtX$XtX
      Xty <- Xty/sdx
    }else{
      sdx <- rep(1,p)
    }

    if(is.null(lambda)){
      Cmax <- ifelse(alpha>0, max(abs(Xty)/alpha), 5)
      Cmax <- ifelse(Cmax <= .Machine$double.eps,1E-5,Cmax)
      lambda <- exp(seq(log(Cmax),log(1E-5),length=nLambda))
      lambda[nLambda] <- 0
    }
    nLambda <- length(lambda)

    beta <- .Call('updatebeta_lambda',as.integer(p),XtX,Xty,as.integer(nLambda),
      as.numeric(lambda),as.numeric(alpha),as.numeric(tol),as.integer(maxIter),verbose)[[1]]

    if(scale) beta <- scale(beta,FALSE,sdx)
    df <- apply(beta,1,function(x)sum(abs(x)>0))

    out <- list(beta=Matrix::Matrix(beta,sparse=TRUE),lambda=lambda,df=df,sdx=sdx)
    return(out)
}
