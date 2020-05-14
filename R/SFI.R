# X = Z = indexK = h2 = subset = lambda = saveAt = name = NULL
# trn =  tst = seq_along(y); alpha = 1; nLambda = 100, method = "CD1"
# mc.cores = 5; tol = 1E-4; maxIter = 500; verbose = TRUE

SFI <- function(y, X = NULL, b = NULL, Z = NULL, K, indexK = NULL,
         h2 = NULL, trn = seq_along(y), tst = seq_along(y), subset = NULL,
         alpha = 1, lambda = NULL, nLambda = 100, method = c("CD1","CD2","LAR"),
         tol = 1E-4, maxIter = 500, maxDF = 100, saveAt = NULL, name = NULL,
         mc.cores = getOption("mc.cores", 2L), verbose = TRUE)
{
  method <- match.arg(method)

  n <- length(y);  nTRN <- length(trn);  nTST <- length(tst)
  if(is.character(K)){
      K <- readBinary(K,indexRow=indexK,indexCol=indexK)
  }
  if(!is.double(K)) K <- apply(K,2,as.double)

  if(is.null(X))   # Design matrix for fixed effects including the intercept
  {
    X <- model.matrix(~1,data=data.frame(rep(1,n)))
  }else{
    if(is.vector(X)){
      X <- stats::model.matrix(~X)
      if(ncol(X)>2)  colnames(X)[-1] <- substr(colnames(X)[-1],2,nchar(colnames(X)[-1]))
    }
  }

  if(!is.null(Z)){
    if(!is.matrix(Z)) stop("Object 'Z' must be a matrix with nrow(Z)=n and ncol(Z)=nrow(K)\n")
    K <- Z%*%K%*%t(Z)
  }

  if(!is.matrix(K) | (length(K) != n^2))
    stop("Product Z %*% K %*% t(Z) must be a squared matrix with number of rows",
     "\n(and columns) equal to the number of elements in 'y'")

  RHS <- K[trn,tst,drop=FALSE]
  K <- K[trn,trn]

  if(is.null(h2))
  {
    # Fit LMM to get variance components and estimate fixed effects as GLS
    fm <- solveMixed(y[trn],X=X[trn,,drop=FALSE],K=K)
    if(fm$convergence){
      varU <- fm$varU
      varE <- fm$varE
      h2 <- varU/(varU + varE)
      b <- fm$b
    }else stop("Convergence was not reached in the 'EMMA' algorithm. \n\t",
           "Please provide a heritability estimate in 'h2' parameter")
  }else{   # Only estimate fixed effects as GLS
    varU <- varE <- NULL
    if(is.null(b)){
      b <- solveMixed(y[trn],X=X[trn,,drop=FALSE],K=K,BLUP=FALSE,h2=h2)$b
    }else{
      if(length(b) != ncol(X)) stop("The length of 'b' must be the same as the number of columns of 'X'\n")
    }
  }
  if(h2 < 0.001) warning("The 'heritability' is too small. Results may be affected",immediate.=TRUE)
  Xb <- drop(X%*%b)

  # Standardizing
  for(i in 1:nTRN)  K[i,i] <- K[i,i] + (1-h2)/h2
  sdx <- sqrt(diag(K))
  scale_cov(K,void=TRUE)   # Equal to K=cov2cor(K) but faster
  RHS <- apply(RHS,2,function(x)x/sdx)

  if(method %in% c("CD1","CD2")){
    if(is.null(lambda)){
      if(method == "CD1"){
          Cmax <- ifelse(alpha > .Machine$double.eps,max(abs(RHS)/alpha),5)
          lambda <- exp(seq(log(Cmax),log(.Machine$double.eps^0.5),length=nLambda))
      }
    }else{
      if(is.matrix(lambda)){
          if(nrow(lambda) != nTST) stop("Object 'lambda' must be a vector or a matrix with nTST rows")
      }
    }
  }

  name <- ifelse(is.null(name),method,name)

  # Split the testing set into subsets. Only the subset provided will be fitted
  if(!is.null(subset)){
     if(!is.numeric(subset) & length(subset) != 2)
      stop("Object 'subset' must contain at least a 2-elements vector")
     partitions <- sort(rep(1:subset[2],ceiling(nTST/subset[2]))[1:nTST])
     index <- which(partitions == subset[1])
     tmp <- paste0(" of ",length(tst))
     tst <- tst[index]
     RHS <- RHS[,index,drop=FALSE]
  }else tmp <- ""

  if(verbose)
    cat(" Fitting SFI model for nTST=",length(tst),tmp," and nTRN=",nTRN," individuals\n",sep="")

  compApply <- function(chunk)
  {
    rhs <- drop(RHS[,chunk])

    if(method == "LAR"){
      fm <- lars2(K,rhs,method="LAR",maxDF=maxDF,scale=FALSE)
      fm$beta <- fm$beta[-1,]
      fm$df <- fm$df[-1]
      fm$lambda <- fm$lambda[-length(fm$lambda)]  
    }else{
      if(is.matrix(lambda)){
        lambda0 <- lambda[chunk,]
      }else lambda0 <- lambda

      fm <- solveEN(K,rhs,scale=FALSE,lambda=lambda0,nLambda=nLambda,
                    alpha=alpha,tol=tol,maxIter=maxIter)
    }

    # Return betas to the original scale by dividing by sdx
    B <- Matrix::Matrix(scale(fm$beta,FALSE,sdx), sparse=TRUE)

    if(verbose){
      cat(1,file=con,append=TRUE)
      utils::setTxtProgressBar(pb, nchar(scan(con,what="character",quiet=TRUE))/length(tst))
    }
    return(list(B=B,lambda=fm$lambda,tst=tst[chunk],df=fm$df))
  }

  if(verbose){
     pb = utils::txtProgressBar(style=3)
     con <- tempfile()
  }
  if(mc.cores == 1L) {
    out = lapply(X=seq_along(tst),FUN=compApply)
  }else{
    out = mclapply(X=seq_along(tst),FUN=compApply,mc.cores=mc.cores)
  }
  if(verbose) {
    close(pb); unlink(con)
  }

  if(sum(tst != unlist(lapply(out,function(x)x$tst)))>0){
    stop("Some sub-processes failed. Something went wrong during the analysis.")
  }

  pMin <- min(unlist(lapply(out,function(x)length(x$lambda))))
  out <- list(method=method, name=name, y=y, Xb=Xb, b=b, varU=varU,
              varE=varE, h2=h2, trn=trn, tst=tst, alpha=alpha,
              df = do.call("rbind",lapply(out,function(x)x$df[1:pMin])),
              lambda = do.call("rbind",lapply(out,function(x)x$lambda[1:pMin])),
              BETA = lapply(out,function(x) x$B[1:pMin, ,drop=FALSE]))
  class(out) <- "SFI"

  # Save outputs if 'saveAt' is not NULL
  if(!is.null(saveAt)){
    if(!is.null(subset)){
       filenames <- paste0(saveAt,c("_B_","_"),subset[1],"_of_",subset[2],c(".bin",".RData"))
    }else  filenames <- paste0(saveAt,c("_B.bin",".RData"))

    if(!file.exists(dirname(filenames[1])))    dir.create(dirname(filenames[1]),recursive = TRUE)
    saveBinary(as.matrix(do.call(rbind,out$BETA)),filename=filenames[1],size=4,verbose=FALSE)
    out$BETA <- filenames[1]
    save(out,file=filenames[2])
  }
  return(out)
}
