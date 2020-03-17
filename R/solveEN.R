
solveEN <- function(P, cov, alpha = 1, lambda = NULL, nLambda = 100,
        scale = TRUE, tol = 1E-5, maxIter = 1000, verbose = FALSE)
{
    p <- length(cov)
    if(length(P) != p^2)
      stop("Incompatible dimensions between 'P' and 'cov'")
    if((sum(dim(P))/2)^2 != length(P))
      stop("Object 'P' must be a squared matrix. The algorithm can not be implemented")

    if(alpha<0 | alpha>1) stop("Parameter 'alpha' must be a number between 0 and 1")

    if(!is.double(cov)) cov <- as.double(cov)
    if(!is.double(P)) P <- apply(P,2,as.double)

    if(scale)
    {
      sdx <- sqrt(diag(P))
      P <- scale_cov(P)  # Equal to P=cov2cor(P) but faster
      cov <- cov/sdx
    }else{
      sdx <- rep(1,p)
    }

    if(is.null(lambda)){
      Cmax <- ifelse(alpha > .Machine$double.eps, max(abs(cov)/alpha), 5)
      lambda <- exp(seq(log(Cmax),log(.Machine$double.eps^0.5),length=nLambda))
    }
    nLambda <- length(lambda)

    beta <- .Call('updatebeta',as.integer(p),P,as.vector(cov),
               as.integer(nLambda),as.numeric(lambda),as.numeric(alpha),
               as.numeric(tol),as.integer(maxIter),verbose)[[1]]

    if(scale) beta <- scale(beta,FALSE,sdx)
    df <- apply(beta,1,function(x)sum(abs(x)>0))

    out <- list(beta=Matrix::Matrix(beta,sparse=TRUE),lambda=lambda,df=df,sdx=sdx)
    class(out) <- "SSI"
    return(out)
}
