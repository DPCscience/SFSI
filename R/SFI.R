
SFI <- function(G,y,h2=0.5,trn=1:length(y),tst=1:length(y),indexG=NULL,subset=NULL,
    kernel=NULL,lambda=NULL,nLambda=100,method=c("CD1","CD2"),alpha=1,name=NULL,
    mc.cores=getOption("mc.cores", 2L),tol=1E-5,maxIter=800,saveAt=NULL,verbose=TRUE)
{
  method <- match.arg(method)

  nTRN <- length(trn);  nTST <- length(tst)
  if(is.character(G)){
      G <- readBinary(G,indexRow=indexG,indexCol=indexG)
  }

  if(!is.matrix(G) | (length(G) != length(y)^2))
    stop("Object 'G' must be a squared matrix with number of rows (and columns) equal to the number of elements in 'y'")

  if(!is.double(G)) G <- apply(G,2,as.double)

  RHS <- G[trn,tst,drop=FALSE]
  G <- G[trn,trn]

  for(i in 1:nTRN)  G[i,i] <- G[i,i] + (1-h2)/h2

  if(!is.null(kernel)){
    if(is.list(kernel) & is.null(kernel$kernel)) stop("Parameter 'kernel' must be a 'list' type object")
    kernel2(G,kernel,void=TRUE)
  }

  # Standardizing
  sdx <- sqrt(diag(G))
  scale_cov(G,void=TRUE)
  RHS <- apply(RHS,2,function(x)x/sdx)

  if(is.null(lambda)){
    if(method == "CD1"){
        Cmax <- ifelse(alpha>0,max(abs(RHS)/alpha),5)
        lambda <- exp(seq(log(Cmax),log(1E-5),length=nLambda))
        lambda[nLambda] <- 0
    }
  }else{
    if(is.matrix(lambda)){
        if(nrow(lambda) != nTST) stop("Object 'lambda' must be a vector or a matrix with nTST rows")
    }
  }

  if(is.null(name)) name <- paste(c(method,substr(kernel$kernel,1,2)),collapse="_")

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

  cat(" Fitting SFI model for nTST=",length(tst),tmp," and nTRN=",length(trn)," individuals\n",sep="")

  compApply <- function(chunk)
  {
    rhs <- drop(RHS[,chunk])
    if(is.matrix(lambda)){ 
      lambda0 <- lambda[chunk,]
    }else lambda0 <- lambda
      
    fm <- solveEN(G,rhs,scale=FALSE,lambda=lambda0,nLambda=nLambda,
                  alpha=alpha,tol=tol,maxIter=maxIter)

    # Returning betas to their original scale by scaling them by their SD
    B <- Matrix::Matrix(scale(fm$beta,FALSE,sdx), sparse=TRUE)

    cat(1,file=con,append=TRUE)
    if(verbose){
      utils::setTxtProgressBar(pb, nchar(scan(con,what="character",quiet=TRUE))/length(tst))
    }
    return(list(B=B,lambda=fm$lambda,tst=tst[chunk],df=fm$df))
  }

  pb = utils::txtProgressBar(style=3)
  con <- tempfile()
  if(mc.cores == 1L) {
    out = lapply(X=seq_along(tst),FUN=compApply)
  }else{
    out = parallel::mclapply(X=seq_along(tst),FUN=compApply,mc.cores=mc.cores)
  }
  close(pb); unlink(con)

  if(sum(tst != unlist(lapply(out,function(x)x$tst)))>0){
    stop("Matching error. Something went wrong during the analysis.")
  }

  out <- list(method=method, kernel=kernel, name=name, y=y, h2=h2,
              trn=trn,tst=tst, alpha=alpha,
              df=do.call("rbind",lapply(out,function(x)x$df)),
              lambda=do.call("rbind",lapply(out,function(x)x$lambda)),
              BETA=lapply(out,function(x) x$B)
            )
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
