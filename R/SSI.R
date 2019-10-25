
SSI <- function(XtX,Xty,kernel=NULL,scale=TRUE,maxDF=NULL,lambda=NULL,nLambda=100,
  method=c("CD","LAR","LAR-LASSO"),alpha=1,name=NULL,tol=1E-5,maxIter=1000,verbose=FALSE)
{
  method <- match.arg(method)

  if(!is.null(dim(Xty)) | !is.double(Xty))  Xty <- as.double(Xty)
  if(!is.double(XtX)) XtX <- apply(XtX,2,as.double)

  if(!is.null(kernel)){
    if(is.list(kernel) & is.null(kernel$kernel)) stop("Parameter 'kernel' must be a 'list' type object")
    XtX <- kernel2(XtX,kernel)
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
