#====================================================================
#====================================================================
cov2dist <- function(XtX,void=FALSE)
{
    if((sum(dim(XtX))/2)^2 != length(XtX)) stop("Object 'XtX' must be a squared matrix")
    if(!is.double(XtX)) XtX <- apply(XtX,2,as.double)
    p <- ncol(XtX)

    if(void){
      tmp <- .Call('crossprod2distance',as.integer(p),XtX)
    }else{
      tmp <- XtX[]
      .Call('crossprod2distance',as.integer(p),tmp)
      tmp
    }
}

#====================================================================
#====================================================================
kernel2 <- function(XtX,kernel=NULL,void=FALSE)
{
    kernel$kernel <- match.arg(kernel$kernel,choices=c("GAUSSIAN","LAPLACIAN","POLYNOMIAL"))

    if((sum(dim(XtX))/2)^2 != length(XtX))
      stop("Object 'XtX' must be a squared matrix. Kernel can not be implemented")
    if(!is.double(XtX)) XtX <- apply(XtX,2,as.double)
    p <- ncol(XtX)

    if(kernel$kernel=="GAUSSIAN"){
        if(is.null(kernel$h)){
          kernel$h <- -2*log(0.5)
          cat("Hyperparameter 'h' was set to ",kernel$h," for GAUSSIAN kernel \n",sep="")
        }else if(kernel$h <= 0) stop("Parameter 'h' for GAUSSIAN kernel must be greater than zero")

        if(void){
          tmp <- .Call('gaussian_kernel',as.integer(p),XtX,as.numeric(kernel$h))
        }else{
          tmp <- XtX[]
          .Call('gaussian_kernel',as.integer(p),tmp,as.numeric(kernel$h))
          tmp
        }

    }else if(kernel$kernel=="LAPLACIAN"){
        if(is.null(kernel$h)){
          kernel$h <- -2*log(0.5)
          cat("Hyperparameter 'h' was set to ",kernel$h," for LAPLACIAN kernel \n",sep="")
        }else if(kernel$h <= 0) stop("Parameter 'h' for LAPLACIAN kernel must be greater than zero")

        if(void){
          tmp <- .Call("laplacian_kernel",as.integer(p),XtX,as.numeric(kernel$h))
        }else{
          tmp <- XtX[]
          .Call("laplacian_kernel",as.integer(p),tmp,as.numeric(kernel$h))
          tmp
        }

    }else if(kernel$kernel=="POLYNOMIAL"){
        if(is.null(kernel$a)){
          kernel$a <- 1
          cat("Hyperparameter 'a' was set to ",kernel$a," for POLYNOMIAL kernel \n",sep="")
        }
        if(is.null(kernel$b)){
          kernel$b <- 2
          cat("Hyperparameter 'b' was set to ",kernel$b," for POLYNOMIAL kernel \n",sep="")
        }else if(round(kernel$b)!=kernel$b) stop("Parameter 'b' for POLYNOMIAL kernel must be an integer")

        if(void){
          tmp <- .Call("polynomial_kernel",as.integer(p),XtX,as.numeric(kernel$a),as.numeric(kernel$b))
        }else{
          tmp <- XtX[]
          .Call("polynomial_kernel",as.integer(p),tmp,as.numeric(kernel$a),as.numeric(kernel$b))
          tmp
        }

    }else{
        stop("Parameter 'kernel' does not fulfill the requirements")
    }
}

#====================================================================
#====================================================================
scale_cov <- function(XtX,void=FALSE)
{
    if((sum(dim(XtX))/2)^2 != length(XtX))
      stop("Object 'XtX' must be a squared matrix. Scaling can not be implemented")
    if(!is.double(XtX)) XtX <- apply(XtX,2,as.double)
    p <- ncol(XtX)

    if(void){
      tmp <- .Call('scaleXtX',as.integer(p),XtX)
    }else{
      tmp <- XtX[]
      .Call('scaleXtX',as.integer(p),tmp)
      tmp
    }
}

#====================================================================
#====================================================================
getIndexCorrelated <- function(X,maxCor=0.8)
{
  COV <- stats::cov(X)
  p <- ncol(COV)
  index <- .Call("getCorrelated",as.integer(p),COV,as.numeric(maxCor))
  out <- NULL
  if(index[[2]]>0) out <- index[[1]][1:index[[2]]]
  out
}

#====================================================================
#====================================================================
collect <- function(prefix="")
{
  filenames <- Sys.glob(paste0(prefix,"_*_of_*.RData"))
  out <- NULL
  if(length(filenames)>0){
      nFiles <- as.numeric(unlist(lapply(strsplit(filenames,"_"),function(x) gsub(".RData","",x[length(x)]))))
      if(length(unique(nFiles))>1)
        stop(" Different subset output files were found for the given prefix='",prefix,
            "'. Remove old files. No output was collected")

      filenames <- paste0(prefix,"_",1:nFiles[1],"_of_",nFiles[1],".RData")
      if(!all(file.exists(filenames))) stop("Some files are missing for the given prefix='",prefix,"'\n")
      nset <- rep(NA,length(filenames))

      for(i in seq_along(filenames))
      {
        load(filenames[i])
        nset[i] <- length(out$tst)
        if(i==1){
          fm <- out
        }else{
          fm$BETA <- c(fm$BETA,out$BETA)
          fm$tst <- c(fm$tst,out$tst)
          fm$df <- rbind(fm$df,out$df)
          fm$lambda <- rbind(fm$lambda,out$lambda)
        }
        cat(" Loaded file: '",filenames[i],"'\n",sep="")
      }
      fm$nset <- nset

  }else stop(" No output files were found for the given prefix='",prefix,"'")

  fm
}

#====================================================================
#====================================================================
backsolvet <- function(r, x, k=ncol(r))
{
  backsolve(r,x,k,transpose=TRUE)
}

#====================================================================
#====================================================================
upDateR <- function(xtx, R = NULL, Xtx, eps = .Machine$double.eps)
{
  norm.xnew <- sqrt(xtx)
  if(is.null(R)) {
    R <- matrix(norm.xnew, 1, 1)
    attr(R, "rank") <- 1
    return(R)
  }
  r <- backsolvet(R, Xtx)
  rpp <- norm.xnew^2 - sum(r^2)
  rank <- attr(R, "rank")	### check if R is machine singular
  if(rpp <= eps)
    rpp <- eps
  else {
    rpp <- sqrt(rpp)
    rank <- rank + 1
  }
  R <- cbind(rbind(R, 0), c(r, rpp))
  attr(R, "rank") <- rank
  R
}

#====================================================================
#====================================================================
downDateR <- function(R, k = p)
{
	p <- dim(R)[1]
	if(p == 1)
		return(NULL)
	R <- deleteCol(R, rep(1, p), k)[[1]][ - p,  , drop = FALSE]
	attr(R, "rank") <- p - 1
	R
}

#====================================================================
#====================================================================
deleteCol <- function(r, z, k = p)
{
	p <- dim(r)[1]
	r <- r[,  - k, drop = FALSE]
	z <- as.matrix(z)
	dz <- dim(z)
	storage.mode(r) <- storage.mode(z) <- "double"
  .Call("delete_col",r,as.integer(p),as.integer(k),z,as.integer(dz[2]))
}

#====================================================================
#====================================================================
saveBinary <- function(X,filename=paste0(tempdir(),"/file.bin"),size=4,verbose=TRUE)
{
  if(!is.matrix(X)) stop("Object 'X' must be a matrix")
  if(storage.mode(X) != "double") storage.mode(X) <- "double"

  unlink(filename)
  out <- .Call('writeBinFile',filename,nrow(X),ncol(X),
          as.integer(size),X)
   if(verbose) cat("Saved file '",filename,"'\n nrows=",nrow(X)," ncols=",nrow(X)," size=",size," bytes\n")
}

#====================================================================
#====================================================================
readBinary <- function(filename=paste0(tempdir(),"/file.bin"),indexRow=NULL,indexCol=NULL,verbose=TRUE)
{
  if(!file.exists(filename)){
     stop("File '",filename,"' does not exist")
  }

  nsetRow <- as.integer(length(indexRow))
  nsetCol <- as.integer(length(indexCol))

  # Read lines
  X <- .Call("readBinFile",filename,nsetRow,nsetCol,
      as.integer(indexRow),as.integer(indexCol))
  n <- X[[1]]; p <- X[[2]]; sizevar <- X[[3]];
  X <- X[[4]]
  if(verbose) cat("Loaded file '",filename,"'\n nrows=",n," ncols=",p," size=",sizevar," bytes\n")

  return(X)
}

#====================================================================
# Labels and breaks for the DF axis
#====================================================================
getSecondAxis <- function(lambda,df,maxLength=6)
{
  loglambda <- -log(lambda)
  labels0 <- sort(unique(round(df)))
  if(min(labels0)<1) labels0[which.min(labels0)] <- 1

  if(stats::IQR(df)>0)
  {
    breaks0 <- stats::predict(stats::smooth.spline(df, loglambda),labels0)$y
  }else breaks0 <- NULL

  index <- 1
  while(any((breaks0-breaks0[max(index)])>1)){
    dd <- breaks0-breaks0[max(index)]
    index <- c(index,which(dd > 1)[1])
  }
  breaks0 <- breaks0[index]
  labels0 <- labels0[index]
  
  if(length(breaks0)>maxLength){
    index <- unique(round(seq(1,length(breaks0),length=maxLength)))
    breaks0 <- breaks0[index]
    labels0 <- labels0[index]
  }
    
  return(list(breaks=breaks0,labels=labels0))
}

#====================================================================
#====================================================================
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
  |====================================================================|
  |    ._______. ._______. ._______. ._______.                         |
  |    | ._____| | ._____| | ._____| |__. .__|                         |
  |    | |_____. | |___.   | |_____.    | |                            |
  |    |_____. | | .___|   |_____. |    | |    Authors:                |
  |    ._____| | | |       ._____| | .__| |__.  Marco Lopez-Cruz       |
  |    |_______| |_|       |_______| |_______|  Gustavo de los Campos  |
  |                                                                    |
  |  Sparse Family and Selection Index. Version 0.1.0 (Oct 31-2019)    |
  |====================================================================|
  ")
}
