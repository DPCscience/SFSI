# X = NULL; Z = NULL; K=G2.0; indexK = NULL; b=fm0$b; h2 = h20; trn.CV = indexTRN; alpha = 1;
# lambda = NULL; nLambda = 100; nCV = 2; nFolds = 5; seed = NULL; method = c("CD1","CD2")[1];
# mc.cores = 4; tol = 1E-4; maxIter = 500; name = NULL; verbose = TRUE

SFI_CV <- function(y, X = NULL, b = NULL, Z = NULL, K, indexK = NULL,
                   h2 = NULL, trn.CV = seq_along(y), alpha = 1, lambda = NULL,
                   nLambda = 100, nCV = 1, nFolds = 5, seed = NULL,
                   method = c("CD1","CD2"), tol = 1E-4, maxIter = 500, name = NULL,
                   lowertri = FALSE, mc.cores = 1, verbose = TRUE)
{
    nFolds <- match.arg(as.character(nFolds),choices=c(2,3,4,5,10,'n'))
    method <- match.arg(method)

    if(is.character(K)){
        K <- readBinary(K,indexRow=indexK,indexCol=indexK)
    }

    if(!is.double(K)) K <- apply(K,2,as.double)

    if(!is.null(Z)){
      if(!is.matrix(Z)){
        stop("Object 'Z' must be a matrix with nrow(Z)=n and ncol(Z)=nrow(K)\n")
      }
      K <- Z%*%K%*%t(Z)
    }

    name <- ifelse(is.null(name),method,name)
    nTRN <- length(trn.CV)
    isLOOCV <- nFolds=='n'
    nFolds <- ifelse(isLOOCV,nTRN,as.numeric(nFolds))
    mc.cores2 <- ifelse(isLOOCV,1,mc.cores)

    compApply <- function(chunk)
    {
      trn <- trn.CV[folds != chunk]
      tst <- trn.CV[folds == chunk]

      fm <- SFI(y, X=X, b=b, K=K, h2=h2, trn=trn, tst=tst,
              alpha=alpha, method=method, lambda=lambda,
              nLambda=nLambda, tol=tol, maxIter=maxIter,
              mc.cores=mc.cores2, lowertri=lowertri, verbose=FALSE)
      fv <- summary(fm)

      if(isLOOCV){
          out <- list(u=as.vector(fitted.SFI(fm)), h2=fm$h2, b=fm$b, tst=tst,
                      df=fm$df, lambda=fm$lambda)
      }else{
          fv <- summary(fm)
          out <- list(h2=fm$h2, b=fm$b, tst=tst ,df=fv$df, lambda=fv$lambda,
                      accuracy=fv$accuracy, MSE=fv$MSE)
      }

      if(verbose){
        cat("Cross-validation:",ifelse(isLOOCV,"LOO",k),". Fold ",chunk," of ",nFolds,"\n",sep="")
      }
      out
    }

    if(is.null(seed)){   # Seeds for randomization
      seeds <- seq((runif(1)+0.05)*1000,.Machine$integer.max,length=nCV)
    }else{
      seeds <- seed
      nCV <- length(seeds)
    }
    stopifnot(nCV == length(seeds))
    nCV <- ifelse(isLOOCV & nCV > 1,1, nCV)   # No need of repeating CV when LOO

    res <- vector("list",nCV)
    for(k in 1:nCV)
    {
      # Creating folds
      folds <- rep(seq(1:nFolds),ceiling(nTRN/nFolds))[1:nTRN]
      if(!isLOOCV){
        set.seed(seeds[k])
        folds <- sample(folds)
      }

      if(mc.cores == 1L | !isLOOCV){
        out = lapply(X=seq(nFolds),FUN=compApply)
      }
      if(mc.cores > 1L & isLOOCV){
        out = parallel::mclapply(X=seq(nFolds),FUN=compApply,mc.cores=mc.cores)
      }

      tmp <- do.call(rbind,split(data.frame(trn.CV,folds),folds))[,1]
      if(sum(tmp != unlist(lapply(out,function(x)x$tst)))>0){
        stop("Some sub-processes failed. Something went wrong during the analysis.")
      }

      # Calculate accuracy
      if(isLOOCV){
        uHat <- do.call(rbind,lapply(out,function(x)x$u))
        accuracy <- suppressWarnings(stats::cor(y[trn.CV],uHat,use="pairwise.complete.obs"))
        MSE <- suppressWarnings(t(apply((y[trn.CV]-uHat)^2,2,sum,na.rm=TRUE)/length(trn.CV)))
        df <- t(apply(do.call("rbind",lapply(out,function(x)x$df)),2,mean))
        lambda0 <- t(apply(do.call("rbind",lapply(out,function(x)x$lambda)),2,mean))
        b0 <- t(apply(do.call("rbind",lapply(out,function(x)x$b)),2,mean))
        h20 <- mean(unlist(lapply(out,function(x)x$h2)))

      }else{
        accuracy <- do.call("rbind",lapply(out,function(x)x$accuracy))
        MSE <- do.call("rbind",lapply(out,function(x)x$MSE))
        df <- do.call("rbind",lapply(out,function(x)x$df))
        lambda0 <- do.call("rbind",lapply(out,function(x)x$lambda))
        b0 <- do.call("rbind",lapply(out,function(x)x$b))
        h20 <- unlist(lapply(out,function(x)x$h2))
      }

      res[[k]] <- list(method=method, name=name, folds=data.frame(trn.CV,fold=folds),
                   b=b0, h2=h20, accuracy=accuracy, MSE=MSE, df=df, lambda=lambda0)
    }
    class(res) <- "SFI_CV"
    return(res)
}
