
SFI_CV <- function(G,y,h2=0.5,trn=1:length(y),nFolds=5,indexG=NULL,kernel=NULL,
    lambda=NULL,nLambda=100,method=c("CD1","CD2","GBLUP"),alpha=1,seed=123,
    mc.cores=getOption("mc.cores", 2L),tol=1E-5,maxIter=800,name=NULL,verbose=TRUE)
{
    nFolds <- match.arg(as.character(nFolds),choices=c(2,3,4,5,10))
    method <- match.arg(method)

    if(is.character(G)){
        G <- readBinary(G,indexRow=indexG,indexCol=indexG)
    }

    if(is.null(name)) name <- paste(c(method,substr(kernel$kernel,1,2)),collapse="_")
    nFolds <- as.numeric(nFolds)
    nTRN <- length(trn)

    # Creating folds
    set.seed(seed)
    folds <- rep(seq(1:nFolds),ceiling(nTRN/nFolds))[1:nTRN]
    folds <- sample(folds)

    out <- vector("list",nFolds)
    for(j in 1:nFolds)
    {
        if(verbose) cat(" Fold ",j," of ",nFolds,"\n",sep="")
        if(method == "GBLUP"){
            fm <- GBLUP(G,y,h2,trn[folds!=j],tst=trn[folds==j],kernel=kernel)
        }else{
            fm <- SFI(G,y,h2,trn[folds!=j],tst=trn[folds==j],kernel=kernel,mc.cores=mc.cores,
                  alpha=alpha,method=method,lambda=lambda,nLambda=nLambda,tol=tol,maxIter=maxIter,verbose=verbose)
        }
        fv <- summary(fm)
        out[[j]] <- list(correlation=fv$correlation, accuracy=fv$accuracy, MSE=fv$MSE,df=fv$df,lambda=fv$lambda)
    }

    out <- list(method=method, kernel=kernel, name=name,y=y, h2=h2,trn=trn,
                alpha=alpha, folds=data.frame(trn,fold=folds),
                correlation=do.call("rbind",lapply(out,function(x)x$correlation)),
                accuracy=do.call("rbind",lapply(out,function(x)x$accuracy)),
                MSE=do.call("rbind",lapply(out,function(x)x$MSE)),
                df=do.call("rbind",lapply(out,function(x)x$df)),
                lambda=do.call("rbind",lapply(out,function(x)x$lambda))
              )
    class(out) <- "SFI_CV"
    return(out)
}
