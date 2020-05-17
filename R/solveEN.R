
# P=tt2; cov=Xty; alpha = 1; lambda = NULL; nLambda = 100; scale = TRUE
# tol = 1E-5; maxIter = 1000; verbose = FALSE; lowertri=FASLSE

solveEN <- function(P, cov, alpha = 1, lambda = NULL,
       nLambda = 100, scale = TRUE, tol = 1E-5, maxIter = 1000,
       verbose = FALSE, lowertri=FALSE)
{
    p <- length(cov)
    if(length(P) != p^2){
      if(lowertri){
        if(length(P) != p*(p+1)/2)
          stop("'P' must be of length p*(p+1)/2 when 'lowertri=TRUE'")
      }else stop("'P' must be of length p*p where p=length(cov)")
    }

    isVectorized <- (length(P) == p*(p+1)/2) & lowertri

    if(alpha<0 | alpha>1) stop("Parameter 'alpha' must be a number between 0 and 1")

    if(!is.double(cov)) cov <- as.double(cov)
    if(!is.double(P)) P <- apply(P,2,as.double)

    if(scale)
    {
      if(!isVectorized){
        sdx <- sqrt(diag(P))
        P <- scale_cov(P)  # Equal to P=cov2cor(P) but faster
        cov <- cov/sdx
      }else{
        warning("scaling is not implemented when 'P' is provided vectorized")
        sdx <- rep(1,p)
      }
    }else{
      sdx <- rep(1,p)
    }

    if(is.null(lambda)){
      Cmax <- ifelse(alpha > .Machine$double.eps, max(abs(cov)/alpha), 5)
      lambda <- exp(seq(log(Cmax),log(.Machine$double.eps^0.5),length=nLambda))
    }
    nLambda <- length(lambda)

    if(lowertri)
    {
      if(!isVectorized){
        P2 <- P[,1]
        for(j in 2:p) P2 <- c(P2, P[j:p,j])
        P <- P2
      }
      beta <- .Call('updatebeta_lowertri',as.integer(p),P,as.vector(cov),
               as.integer(nLambda),as.numeric(lambda),as.numeric(alpha),
               as.numeric(tol),as.integer(maxIter),verbose)[[1]]

    }else{
      beta <- .Call('updatebeta',as.integer(p),P,as.vector(cov),
               as.integer(nLambda),as.numeric(lambda),as.numeric(alpha),
               as.numeric(tol),as.integer(maxIter),verbose)[[1]]
    }

    if(scale) beta <- scale(beta,FALSE,sdx)
    df <- apply(beta,1,function(x)sum(abs(x)>0))

    out <- list(beta=Matrix::Matrix(beta,sparse=TRUE),lambda=lambda,df=df,sdx=sdx)
    class(out) <- "SSI"
    return(out)
}
