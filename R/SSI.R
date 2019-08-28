#' Sparse Selection Index
#'
#' Computes the entire Elastic-Net solution for the regression coefficients of a penalized regression simultaneously for all
#' values of the penalization parameter via either the Coordinate Descent (Friedman, 2007) or Least Angle Regression (Efron, 2004) algorithms.
#' Uses either 'solveEN' or 'lars2' functions
#'
#' The regression coefficients \eqn{beta=(beta_1,...,beta_p)} are estimated as function of the 'variance' matrix among
#' predictors (XtX) and the 'covariance' vector between response and predictors (Xty) by optimizing the function
#' \deqn{-Xty' beta + 0.5 beta' XtX beta + lambda J(beta)}
#' where \eqn{lambda} is the penalization parameter and \eqn{J(beta)} is a penalty function given by
#' \deqn{0.5(1-alpha)||beta||_2^2 + alpha||beta||_1}
#' @return  List object containing the elements:
#' \itemize{
#'   \item beta: vector of regression coefficients.
#'   \item alpha: value for the elastic-net weights used.
#'   \item lambda: sequence of values of lambda used.
#'   \item df: degrees of freedom, number of non-zero predictors at each solution.
#'   \item sdx: vector of standard deviation of predictors.
#'   \item kernel: transformation applied to the elements of XtX.
#' }
#' The returned object is of the class 'SSI' for which methods 'predict', 'plot' and 'summary' exist
#' @param XtX Variance-covariance matrix among predictors
#' @param Xty Covariance vector between response variable and predictors
#' @param kernel Kernel transformation applied to XtX. List consisting on one of:
#' \itemize{
#'   \item list(kernel='GAUSSIAN',h). If \eqn{h} is not provided the value of \eqn{h=-2*log(0.5)} is used.
#'   \item list(kernel='LAPLACIAN',h). If \eqn{h} is not provided the value of \eqn{h=-2*log(0.5)} is used.
#'   \item list(kernel='POLYNOMIAL',a,b). The values of \eqn{a=1} and \eqn{b=2} are used when they are not provided.
#' }
#' Default kernel=NULL (no kernel)
#' @param method One of:
#' \itemize{
#'  \item 'CD': Coordinate Descent algorithm that computes the coefficients for a provided grid of lambdas.
#'  \item 'LAR': Least Angle Regression algorithm that computes the entire sequence of all coefficients. Values of lambdas are calculated at each step.
#'  \item 'LAR-LASSO': Similar to 'LAR' but solutions when a predictor leaves the solution are also returned.
#' }
#' @param maxDF Maximum (average across individuals) number of predictors in the last solution (when method='LAR' or 'LAR-LASSO').
#' Default \eqn{maxDF=NULL} will calculate solutions including 1,2,...,nTRN predictors
#' @param lambda Penalization parameter sequence vector used for the Coordinate Descent algorithm.
#' Default is \eqn{lambda=NULL}, in this case a decreasing grid of
#' n='nLambda' lambdas will be generated starting from a maximum equal to \eqn{max(abs(Xty)/alpha)} to a minumum equal to zero.
#' If \eqn{alpha=0} the grid is generated starting from a maximum equal to 5
#' @param nLambda Number of lambdas generated when \eqn{lambda=NULL}
#' @param alpha Numeric between 0 and 1 indicating the weights for LASSO (alpha) and Ridge-Regression (1-alpha)
#' @param scale TRUE or FALSE to whether scaling each entry of XtX and Xty
#' by the SD of the corresponding predictor taken from the diagonal of XtX
#' @param tol Maximum error between two consecutive solutions of the iterative algorithm to declare convergence
#' @param maxIter Maximum number of iterations to run at each lambda step before convergence is reached
#' @param name Name given to the output for tagging purposes. Default \eqn{name=NULL} will give the name of the method used
#' @param verbose TRUE or FALSE to whether printing each step
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate variables
#' n = 500; p=200;  rho=0.65
#' X = matrix(rnorm(n*p),ncol=p)
#' eta = scale(X%*%rnorm(p))  # signal
#' e =  rnorm(n)              # error
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
#' fm = SSI(P,rhs)                    # or
#' fm = SSI(P,rhs,method="LAR")       # or
#' fm = SSI(P,rhs,method="LAR-LASSO")
#'
#' # Regression coefficients
#' plot(fm)  # Path plot along lambda
#' beta = as.matrix(fm$beta)
#'
#' # Predicted values in training and testing set
#' yHat_TRN =  X[trn,] %*% t(beta)   # or
#' yHat_TRN =  predict(fm,X[trn,])
#'
#' yHat_TST =  X[tst,] %*% t(beta)   # or
#' yHat_TST =  predict(fm,X[tst,])
#'
#' par(mfrow=c(1,2))
#' plot(fm$df,cor(y[trn],yHat_TRN)[1,],main="Training set")
#' plot(fm$df,cor(y[tst],yHat_TST)[1,],main="Testing set")
#' @export
#' @references
#' \itemize{
#' \item \insertRef{Efron2004}{SFSI}
#' \item \insertRef{Friedman2007}{SFSI}
#' \item \insertRef{Hoerl1970}{SFSI}
#' \item \insertRef{Tibshirani1996}{SFSI}
#' \item \insertRef{VanRaden2008}{SFSI}
#' \item \insertRef{Zou2005}{SFSI}
#' }
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#' @keywords SSI

SSI <- function(XtX,Xty,kernel=NULL,scale=TRUE,maxDF=NULL,lambda=NULL,nLambda=100,
  method=c("CD","LAR","LAR-LASSO"),alpha=1,name=NULL,tol=1E-5,maxIter=1000,verbose=FALSE)
{
  method <- match.arg(method)

  if(!is.null(dim(Xty)) | !is.double(Xty))  Xty <- as.double(Xty)
  if(!is.double(XtX)) XtX <- apply(XtX,2,as.double)

  if(!is.null(kernel)){
    if(is.list(kernel) & is.null(kernel$kernel)) stop("Parameter 'kernel' must be a 'list' type object")
    XtX <- kernel2(XtX,kernel)
    kernel <- XtX$kernel
    XtX <- XtX$K
  }

  if(!is.null(lambda) | alpha!=1){
    if(method != "CD") cat("Method is changed to 'CD' when lambda is provided and/or alpha<1\n")
    method <- "CD"
  }

  if(is.null(name)) name <- paste(c(method,substr(kernel$kernel,1,2)),collapse="_")

  if(method %in% c("LAR","LAR-LASSO")){
    fm <- lars2(XtX,Xty,method=method,maxDF=maxDF,scale=scale,verbose=verbose)
  }
  if(method == "CD"){
    fm <- solveEN(XtX,Xty,nLambda=nLambda,lambda=lambda,alpha=alpha,
        scale=scale,tol=tol,maxIter=maxIter,verbose=verbose)
  }

  out <- list(alpha=alpha, lambda=fm$lambda, df=fm$df, beta=fm$beta,
    sdx=fm$sdx, method=method, kernel=kernel, name=name)
  class(out) <- "SSI"
  return(out)
}
