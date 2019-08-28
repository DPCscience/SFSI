#====================================================================
#' cov2dist function
#'
#' Computes the pairwise squared Euclidean distance \eqn{d(x,y)} from the covariance between variables 'x' and 'y'
#' given by the cross-product \eqn{x'y}
#' \deqn{d^2(x,y)=(x-y)'(x-y)=x'x + y'y - 2x'y}
#' If the vectors 'x' and 'y' are centered then x'y/(n-1)=cov(x,y) is the sample covariance, where \eqn{n} is the length of the vectors. Then
#' \deqn{d^2(x,y)/(n-1)=var(x) + var(y) - 2cov(x,y)}
#' @return  A squared matrix D containing the squared Euclidean distances
#' @param XtX Cross-product (\eqn{X'X}) of a matrix \eqn{X} containing \eqn{n} observations and \eqn{p} (centered) variables.
#' If a variance-covariance matrix is provided as \eqn{X'X/(n-1)} then the output is \eqn{D/(n-1)}.
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate matrix
#' n = 100; p=10
#' X = scale(matrix(rnorm(n*p),ncol=p))
#'
#' # Distance matrix from a cross-product
#' COV = crossprod(X)   # Cross-product X'X
#' cov2dist(COV)  # equal to: as.matrix(dist(t(X)))^2
#'
#' # Distance matrix from a variance-covariance matrix
#' COV = var(X)   # Variance matrix of X
#' (n-1)*cov2dist(COV)  # equal to: as.matrix(dist(t(X)))^2
#' @export
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#====================================================================
cov2dist <- function(XtX)
{
    if((sum(dim(XtX))/2)^2 != length(XtX)) stop("Object 'XtX' must be a squared matrix")
    if(!is.double(XtX)) XtX <- apply(XtX,2,as.double)
    p <- ncol(XtX)

    XtX <- .Call('crossprod2distance',as.integer(p),XtX)[[1]]
    return(XtX)
}

#====================================================================
#' kernel2 function
#'
#' Computes kernel transformation \eqn{K(x,y)} from the covariance between variables 'x' and 'y'
#' given by the cross-product \eqn{x'y}
#' \enumerate{
#'   \item Gaussian kernel. Bandwidth parameter \eqn{h}:
#' \deqn{K(x,y)=exp{-h d^2(x,y)}=exp{-h [x'x + y'y - 2x'y]}}
#'   \item Laplacian kernel. Bandwidth parameter \eqn{h}:
#' \deqn{K(x,y)=exp{-h d(x,y)}=exp{-h sqrt[x'x + y'y - 2x'y]}}
#'   \item Polynomial kernel. Parameters \eqn{a>0} and integer \eqn{b}:
#' \deqn{K(x,y)=[a(x'y) + 1]^b}
#' }
#' @return  A squared matrix K containing the kernel transformation
#' @param XtX Cross-product (\eqn{X'X}) of a matrix \eqn{X} containing \eqn{n} observations and \eqn{p} (centered) variables.
#' @param kernel List consisting on one of:
#' \itemize{
#'   \item list(kernel='GAUSSIAN',h). If \eqn{h} is not provided the value of \eqn{h=-2*log(0.5)} is used.
#'   \item list(kernel='LAPLACIAN',h). If \eqn{h} is not provided the value of \eqn{h=-2*log(0.5)} is used.
#'   \item list(kernel='POLYNOMIAL',a,b). The values of \eqn{a=1} and \eqn{b=2} are used when they are not provided.
#' }
#' Default kernel=NULL (no kernel)
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate matrix
#' n = 1000; p=10
#' X = scale(matrix(rnorm(n*p),ncol=p))
#' COV = crossprod(X)/(n-1)   # Cross-product X'X
#'
#' # Gaussian kernel (h=0.2)
#' h = 0.2
#' K = kernel2(COV,kernel=list(kernel="GAUSSIAN",h=h))$K
#' # it must equal (but faster) to:
#' K2 = exp(-h*as.matrix(dist(t(X)))^2/(n-1))  # or
#' K2 = exp(-h*cov2dist(COV))
#' K;K2
#'
#' # Laplacian kernel (h=0.2)
#' K = kernel2(COV,kernel=list(kernel="LAPLACIAN",h=h))$K
#' # it must equal (but faster) to:
#' K2 = exp(-h*as.matrix(dist(t(X)))/sqrt(n-1))  # or
#' K2 = exp(-h*sqrt(cov2dist(COV)))
#' K;K2
#'
#' # Polynomial kernel (a=1.5, b=2)
#' a = 1.5; b = 2
#' K = kernel2(COV,kernel=list(kernel="POLYNOMIAL",a=a,b=b))$K
#' # it must equal to:
#' K2 = (a*COV + 1)^b
#' K;K2
#' @export
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#====================================================================
kernel2 <- function(XtX,kernel=NULL)
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
        }else if(kernel$h <= 0) stop("Parameter 'h' for gaussian kernel must be greater than zero")

        XtX <- .Call('gaussian_kernel',as.integer(p),XtX,as.numeric(kernel$h))[[1]]

    }else if(kernel$kernel=="LAPLACIAN"){
        if(is.null(kernel$h)){
          kernel$h <- -2*log(0.5)
          cat("Hyperparameter 'h' was set to ",kernel$h," for LAPLACIAN kernel \n",sep="")
        }else if(kernel$h <= 0) stop("Parameter 'h' for laplacian kernel must be greater than zero")

        XtX <- .Call("laplacian_kernel",as.integer(p),XtX,as.numeric(kernel$h))[[1]]

    }else if(kernel$kernel=="POLYNOMIAL"){
        if(is.null(kernel$a)){
          kernel$a <- 1
          cat("Hyperparameter 'a' was set to ",kernel$a," for POLYNOMIAL kernel \n",sep="")
        }
        if(is.null(kernel$b)){
          kernel$b <- 2
          cat("Hyperparameter 'b' was set to ",kernel$b," for polynomial kernel \n",sep="")
        }else if(round(kernel$b)!=kernel$b) stop("Parameter 'b' for polynomial kernel must be an integer")

        XtX <- .Call("polynomial_kernel",as.integer(p),XtX,as.numeric(kernel$a),as.numeric(kernel$b))[[1]]
    }else{
        stop("Parameter 'kernel' does not fulfill the requirements")
    }
    return(list(K=XtX,kernel=kernel))
}

#====================================================================
#' scale_crossprod function
#'
#' Scales the cross-product matrix \eqn{X'X} by dividing each element by its corresponding standard deviations.
#' The entry corresponding to the cross-product \eqn{x'y} of the (centered) variables 'x' and 'y' is divided by the product of the standard deviations
#' \eqn{sqrt(x'x)} and \eqn{sqrt(y'y)}. The scaled matrix will have all diagonal entries equal to 1.
#' @return  A list object containing the elements:
#' \itemize{
#'   \item 'XtX': squared matrix containing the scaled cross-product matrix \eqn{X'X}
#'   \item 'sdX': vector containing the standard deviation.
#' }
#' @param XtX Cross-product (\eqn{X'X}) of a matrix \eqn{X} containing \eqn{n} observations and \eqn{p} (centered) variables.
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate matrix
#' n = 100; p=10
#' X = matrix(rnorm(n*p),ncol=p)
#' COV = var(X)   # Cross-product X'X
#'
#' COV2 = scale_crossprod(COV)
#' COV2$XtX   # equal to: var(scale(X))
#' COV2$sdX   # equal to: apply(X,2,sd)
#' @export
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#====================================================================
scale_crossprod <- function(XtX)
{
    if((sum(dim(XtX))/2)^2 != length(XtX))
      stop("Object 'XtX' must be a squared matrix. Scaling can not be implemented")
    if(!is.double(XtX)) XtX <- apply(XtX,2,as.double)
    p <- ncol(XtX)
    XtX <- .Call('scaleXtX',as.integer(p),XtX)
    return(list(XtX=XtX[[1]],sdX=XtX[[2]]))
}

#====================================================================
#====================================================================
getIndexCorrelated <- function(X,maxCor=0.8)
{
  COV <- cov(X)
  p <- ncol(COV)
  index <- .Call("getCorrelated",as.integer(p),COV,as.numeric(maxCor))
  out <- NULL
  if(index[[2]]>0) out <- index[[1]][1:index[[2]]]
  out
}

#====================================================================
#' getDistantSet function
#'
#' Selects a set of \eqn{n} individuals that are maximally distant. The individuals are selected as the closest ones
#' to the centroid of each cluster that are generated using k-means clustering using the principal components of a relationship matrix.
#' @return  A vector with the index of the maximally distant individuals
#' @param G Relationship matrix among individuals
#' @param PC Principal components from the relationship matrix 'G'. When 'PC' is provided the 'G' is not needed
#' @param n Number of individuals to select
#' @param kernel Kernel transformation applied to 'G'. List consisting on one of:
#' \itemize{
#'   \item list(kernel='GAUSSIAN',h). If \eqn{h} is not provided the value of \eqn{h=-2*log(0.5)} is used.
#'   \item list(kernel='LAPLACIAN',h). If \eqn{h} is not provided the value of \eqn{h=-2*log(0.5)} is used.
#'   \item list(kernel='POLYNOMIAL',a,b). The values of \eqn{a=1} and \eqn{b=2} are used when they are not provided.
#' }
#' Default kernel=NULL (no kernel)
#' @export
#====================================================================
getDistantSet <- function(G=NULL,PC=NULL,n,kernel=NULL)
{
  if(!is.null(kernel) & !is.null(G)){
    if(is.null(PC)){
      if(is.list(kernel) & is.null(kernel$kernel)) stop("Parameter 'kernel' must be a 'list' type object")
      G <- kernel2(G,kernel)
      kernel <- G$kernel
      G <- G$K
    }else{
      cat("Kernel is only applied when no 'PC' is provided\n")
    }
  }
  if(is.null(PC)){
    PC <-  eigs_sym(G, 2)$vectors
  }else{
      if(!is.matrix(PC)) stop("'PC' must be a matrix")
  }

  # Get nTST means
  K <- kmeans(PC, centers = n, nstart = 100)

  # Get the closest point to centroid
  set <- c()
  for(j in 1:n)
  {
    index <- which(K$cluster==j)
    tmp <- PC[index,]
    c0 <- K$centers[rep(j,nrow(tmp)),]
    d0 <- apply((c0 - tmp)^2,1,sum)
    set <- c(set,index[which.min(d0)])
  }
  set <- sort(set)
  set
}

#====================================================================
#' collect function
#'
#' Collects all outputs saved when 'saveAt' parameter was provided to run the SFI when this was splited using 'subset' parameter.
#' @return  An object of the class 'SFI'
#' @param prefix Prefix that was added to the output files name, this may include a path
#' @examples
#' require(SFSI)
#' # Read data from BGLR package
#' data(wheat,package="BGLR")
#' X = scale(wheat.X)
#' G = tcrossprod(X)/ncol(X)    # Genomic relationship matrix
#' y = wheat.Y[,1]              # Response variable
#' h2 = 0.5                     # Heritability (rough estimate)
#'
#' # Training and testing sets
#' set.seed(1234)
#' n = length(y)
#' pTST = 0.3      # percentage to predict
#' tst = sample(1:n,floor(pTST*n))
#' trn = (1:n)[-tst]
#'
#' prefix <- "testSFI"
#' # Run the analysis into 5 subsets and save them at a given prefix
#' fm <- SFI(G,y,h2,trn,tst,subset=c(1,5),saveAt=prefix,nCores=5)
#' fm <- SFI(G,y,h2,trn,tst,subset=c(2,5),saveAt=prefix,nCores=5)
#' fm <- SFI(G,y,h2,trn,tst,subset=c(3,5),saveAt=prefix,nCores=5)
#' fm <- SFI(G,y,h2,trn,tst,subset=c(4,5),saveAt=prefix,nCores=5)
#' fm <- SFI(G,y,h2,trn,tst,subset=c(5,5),saveAt=prefix,nCores=5)
#'
#' # Collect all results after completion
#' fm <- collect(prefix)
#' plot(fm)
#' unlink(paste0(prefix,"*.RData"))   # Remove files
#' unlink(paste0(prefix,"*.bin"))     # Remove files
#' @export
#' @references
#' \itemize{
#' \item \insertRef{Perez2014}{SFSI}
#' }
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#====================================================================
collect <- function(prefix="")
{
  filenames <- Sys.glob(paste0(prefix,"_*_of_*.RData"))
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
        nset[i] <- length(out$testing)
        if(i==1){
          fm <- out
        }else{
          fm$BETA <- c(fm$BETA,out$BETA)
          fm$testing <- c(fm$testing,out$testing)
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
#' saveBinary function
#'
#' Save a fortran-formatted binary file which will contain information about the number of rows, number of columns and the type of precision (single or double).
#' @return  NULL
#' @param filename Name of the binary file to read
#' @param size Size of a real variable in bytes. size=4 (single precision), size=8 (double precision)
#' @param verbose TRUE or FALSE to whether printing file information
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate matrix
#' n = 10; p = 10
#' X = matrix(rnorm(n*p),ncol=p)
#' X
#'
#' # Save matrix
#' saveBinary(X,filename="testMatrix1.bin",size=4)  # as single-precision
#' saveBinary(X,filename="testMatrix2.bin",size=8)  # as double-precision
#'
#' # Read the single-precision matrix
#' X2 = readBinary("testMatrix1.bin")
#' X2
#' sum(abs(X-X2))   # Note the loss of precision
#'
#' # Read the double-precision matrix
#' X2 = readBinary("testMatrix2.bin")
#' X2
#' sum(abs(X-X2))   # No loss of precision
#' @export
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#====================================================================
saveBinary <- function(X,filename=paste0(tempdir(),"/file.bin"),size=4,verbose=TRUE)
{
  if(!is.matrix(X)) stop("Object 'X' must be a matrix")
  if(storage.mode(X) != "double") storage.mode(X) <- "double"

  unlink(filename)
  int_filename <- utf8ToInt(filename)
  if(length(int_filename)>100) stop("'Filename' must be no more than 100 characters long")

  out <- .Fortran('writeBinFile2',int_filename,length(int_filename),nrow(X),ncol(X),
          as.integer(size),X,integer(1))[[7]]

   if(verbose) cat("Saved file '",filename,"'\n nrows=",nrow(X)," ncols=",nrow(X)," size=",size," bytes\n",sep="")
}

#====================================================================
#' readBinary function
#'
#' Read a fortran-formatted binary file which contains information about the number of rows, number of columns and the type of precision (single or double).
#' @return  A numeric matrix whose dimensions and precision are specified in the file
#' @param filename Name of the binary file to read
#' @param indexRow The index of the rows to be read from the file (1,2,...). Default \eqn{indexRow=NULL} will read all the rows
#' @param indexCol The index of the columns to be read from the file (1,2,...). Default \eqn{indexCol=NULL} will read all the columns
#' @param verbose TRUE or FALSE to whether printing file information
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate matrix
#' n = 10; p = 10
#' X = matrix(rnorm(n*p),ncol=p)
#' X
#'
#' # Save matrix
#' saveBinary(X,filename="testMatrix1.bin",size=4)  # as single-precision
#' saveBinary(X,filename="testMatrix2.bin",size=8)  # as double-precision
#'
#' # Read the single-precision matrix
#' X2 = readBinary("testMatrix1.bin")
#' X2
#' sum(abs(X-X2))   # Note the loss of precision
#'
#' # Read the double-precision matrix
#' X2 = readBinary("testMatrix2.bin")
#' X2
#' sum(abs(X-X2))   # No loss of precision
#'
#' # Read specific rows and columns
#' indexRow = c(2,4,5,8,10)
#' indexCol = c(1,2,6,7,10)
#' X2 = readBinary("testMatrix2.bin",indexRow=indexRow,indexCol=indexCol)
#' X2
#' sum(abs(X[indexRow,indexCol]-X2))
#' @export
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#====================================================================
readBinary <- function(filename=paste0(tempdir(),"/file.bin"),indexRow=NULL,indexCol=NULL,verbose=TRUE)
{
  if(!file.exists(filename)){
     stop("File '",filename,"' does not exist")
  }

  nsetRow <- as.integer(length(indexRow))
  nsetCol <- as.integer(length(indexCol))
  int_filename <- utf8ToInt(filename)
  if(length(int_filename)>100) stop("'Filename' must be no more than 100 characters long")

  # Read the file info
  f <- file(filename,open="rb")
  nrows <- readBin(f,n=1,what=integer())
  ncols <- readBin(f,n=1,what=integer())
  sizevar <- readBin(f,n=1,what=integer())
  close(f)
  n <- ifelse(nsetRow>0,nsetRow,nrows)
  p <- ifelse(nsetCol>0,nsetCol,ncols)

  # Read lines
  out <- matrix(0,n,p)
  storage.mode(out) <- "numeric"
  out <- .Fortran("readBinFile2",int_filename,length(int_filename),nsetRow,nsetCol,
      as.integer(indexRow),as.integer(indexCol),integer(1),ncols,sizevar,n,p,out)[[12]]

  if(verbose) cat("Loaded file '",filename,"'\n nrows=",n," ncols=",p," size=",sizevar," bytes\n",sep="")

  return(out)
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
  |  Sparse Family and Selection Indices. Version 1.0.0                |
  |====================================================================|
  ")
}