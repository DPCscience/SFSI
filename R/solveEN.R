
solveEN <- function(XtX, Xty, alpha = 1, lambda = NULL, nLambda = 100,
        scale = TRUE, tol = 1E-5, maxIter = 1000, verbose = FALSE)
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
      sdx <- sqrt(diag(XtX))
      XtX <- scale_cov(XtX)  # Equal to XtX=cov2cor(XtX) but faster
      Xty <- Xty/sdx
    }else{
      sdx <- rep(1,p)
    }

    if(is.null(lambda)){
      Cmax <- ifelse(alpha > .Machine$double.eps, max(abs(Xty)/alpha), 5)
      lambda <- exp(seq(log(Cmax),log(.Machine$double.eps^0.5),length=nLambda))
    }
    nLambda <- length(lambda)

    beta <- .Call('updatebeta',as.integer(p),XtX,as.vector(Xty),
               as.integer(nLambda),as.numeric(lambda),as.numeric(alpha),
               as.numeric(tol),as.integer(maxIter),verbose)[[1]]

    if(scale) beta <- scale(beta,FALSE,sdx)
    df <- apply(beta,1,function(x)sum(abs(x)>0))

    out <- list(beta=Matrix::Matrix(beta,sparse=TRUE),lambda=lambda,df=df,sdx=sdx)
    class(out) <- "SSI"
    return(out)
}
