#====================================================================
#' cov2dist function
#'
#' Computes a squared Euclidean distance matrix from a covariance matrix among \eqn{p} variables. The pairwise
#' distance \eqn{d(x,y)} between vectors \eqn{x=(x_1,...,x_n)'} and \eqn{y=(y_1,...,y_n)'}
#' is obtained from their cross-product, \eqn{x'y = \sum{x_iy_i}}, as
#' \deqn{d^2(x,y) = (x-y)'(x-y) = x'x + y'y - 2x'y}
#' Note that if the variables are centered then the cross-product is proportional to the covariance,
#' this is \eqn{x'y = cov(x,y)} up-to a constant
#' @return  A squared matrix \eqn{D} containing the squared Euclidean distances
#' @param XtX Cross-product (\eqn{X'X}) of a matrix \eqn{X} with \eqn{p} (centered) variables
#' @param void \code{TRUE} or \code{FALSE} to whether return or not return the output.
#' When \code{FALSE} no result is displayed but the input is modified. Default \code{void=FALSE}
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate matrix
#' n = 100; p=10
#' X = scale(matrix(rnorm(n*p),ncol=p))
#'
#' # Distance matrix from a cross-product
#' COV = crossprod(X)   # Cross-product X'X
#' cov2dist(COV)
#' # it must equal (but faster) to:
#' as.matrix(dist(t(X)))^2
#'
#' # Distance matrix from a variance-covariance matrix
#' COV = cov(X)   # Variance matrix of X
#' (n-1)*cov2dist(COV)
#' # it must equal (but faster) to:
#' as.matrix(dist(t(X)))^2
#'
#' # Using void=TRUE
#' cov2dist(COV,void=TRUE)
#' (n-1)*COV   # notice that COV was modified
#' @export
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
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
#' kernel2 function
#'
#' Applies a kernel transformation to a covariance matrix among \eqn{p} predictors. The pairwise transformation
#' \eqn{K(x,y)} for vectors \eqn{x=(x_1,...,x_n)'} and \eqn{y=(y_1,...,y_n)'} is calculated from their cross-product,
#' \eqn{x'y = \sum{x_iy_i}}, as
#' \enumerate{
#'   \item Gaussian kernel. Bandwidth parameter \eqn{h}:
#' \deqn{K(x,y)=exp\{-h d^2(x,y)\}=exp\{-h [x'x + y'y - 2x'y]\}}
#'   \item Laplacian kernel. Bandwidth parameter \eqn{h}:
#' \deqn{K(x,y)=exp\{-h d(x,y)\}=exp\{-h \sqrt{x'x + y'y - 2x'y]}\}}
#'   \item Polynomial kernel. Parameters \eqn{a>0} and integer \eqn{b}:
#' \deqn{K(x,y)=[a(x'y) + 1]^b}
#' }
#' @return  A squared matrix \eqn{K} containing the kernel transformation
#' @param XtX Cross-product (\eqn{X'X}) of a matrix \eqn{X} with \eqn{p} (centered) variables.
#' @param kernel List consisting on one of:
#' \itemize{
#'   \item \code{list(kernel='GAUSSIAN',h)}. If \code{h} is not provided the value of \code{h=-2*log(0.5)} is used.
#'   \item \code{list(kernel='LAPLACIAN',h)}. If \code{h} is not provided the value of \code{h=-2*log(0.5)} is used.
#'   \item \code{list(kernel='POLYNOMIAL',a,b)}. The values of \code{a=1} and \code{b=2} are used when they are not provided.
#' }
#' Default \code{kernel=NULL} (no kernel)
#' @param void \code{TRUE} or \code{FALSE} to whether return or not return the output.
#' When \code{FALSE} no result is displayed but the input is modified. Default \code{void=FALSE}
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate matrix
#' n = 1000; p=10
#' X = scale(matrix(rnorm(n*p),ncol=p))
#' COV = cov(X)
#'
#' # Gaussian kernel (h=0.2)
#' h = 0.2
#' kernel2(COV,kernel=list(kernel="GAUSSIAN",h=h))
#' # it must be equal (but faster) to:
#' exp(-h*as.matrix(dist(t(X)))^2/(n-1))  # or
#' exp(-h*cov2dist(COV))
#'
#' # Laplacian kernel (h=0.2)
#' kernel2(COV,kernel=list(kernel="LAPLACIAN",h=h))
#' # it must be equal (but faster) to:
#' exp(-h*as.matrix(dist(t(X)))/sqrt(n-1))  # or
#' exp(-h*sqrt(cov2dist(COV)))
#'
#' # Polynomial kernel (a=1.5, b=2)
#' a = 1.5; b = 2
#' kernel2(COV,kernel=list(kernel="POLYNOMIAL",a=a,b=b))
#' # it must be equal to:
#' (a*COV + 1)^b
#'
#' # Using void=TRUE
#' kernel2(COV,kernel=list(kernel="POLYNOMIAL",a=a,b=b),void=TRUE)
#' COV   # notice that COV was modified
#' @export
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
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
#' scale_crossprod function
#'
#' Recalculate a variance-covariance matrix among \eqn{p} variables to its equivalent for the variables scaled to have unit variance.
#' The recalculated matrix will contain as off-diagonal entries the covariance between the scaled variables
#' \eqn{x*=x/\sigma_x} and \eqn{y*=y/\sigma_y}, formed by dividing original variables \eqn{x} and \eqn{y} by their
#' standard deviation \eqn{\sigma_x} and \eqn{\sigma_y}, given by
#' \deqn{cov(x*,y*)=cov(x,y)/(\sigma_x\sigma_y)}
#' while in the diagonal the variance will be \eqn{var(x*)=var(x)/\sigma^2_x=1}.
#' @return  A squared matrix with the recalculated variances and covariances
#' @param XtX Variance-covariance matrix among \eqn{p} variables
#' @param void \code{TRUE} or \code{FALSE} to whether return or not return the output.
#' When \code{FALSE} no result is displayed but the input is modified. Default \code{void=FALSE}
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Simulate matrix
#' n = 100; p=10
#' X = matrix(rnorm(n*p),ncol=p)
#' COV = cov(X)   # Cross-product X'X
#'
#' scale_crossprod(COV)
#' # it must be equal to:
#' cov(scale(X))    # or
#' SD=apply(X,2,sd); cov(scale(X,center=FALSE,scale=SD))
#'
#' # Using void=TRUE
#' scale_crossprod(COV,void=TRUE)
#' COV   # notice that COV was modified
#' @export
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#====================================================================
scale_crossprod <- function(XtX,void=FALSE)
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
#' collect function
#'
#' Collects all outputs saved at the provided \code{saveAt} parameter from the SFI analysis when testing data was splited
#' according to \code{subset} parameter.
#' @return  An object of the class 'SFI' for which methods \code{fitted}, \code{predict}, \code{plot} and \code{summary} exist
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
#' fm <- SFI(G,y,h2,trn,tst,subset=c(1,5),saveAt=prefix,mc.cores=5)
#' fm <- SFI(G,y,h2,trn,tst,subset=c(2,5),saveAt=prefix,mc.cores=5)
#' fm <- SFI(G,y,h2,trn,tst,subset=c(3,5),saveAt=prefix,mc.cores=5)
#' fm <- SFI(G,y,h2,trn,tst,subset=c(4,5),saveAt=prefix,mc.cores=5)
#' fm <- SFI(G,y,h2,trn,tst,subset=c(5,5),saveAt=prefix,mc.cores=5)
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
#' Save a fortran-formatted binary file at a defined precision (single or double).
#' @return  NULL
#' @param X Numeric matrix to save
#' @param filename Name that will be given to the binary file
#' @param size Size of a real variable in bytes (\code{size=4} for single precision and \code{size=8} for double precision)
#' that matrix \code{X} will occupy
#' @param verbose \code{TRUE} or \code{FALSE} to whether printing file information
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
  out <- .Call('writeBinFile',filename,nrow(X),ncol(X),
          as.integer(size),X)
   if(verbose) cat("Saved file '",filename,"'\n nrows=",nrow(X)," ncols=",nrow(X)," size=",size," bytes\n")
}

#====================================================================
#' readBinary function
#'
#' Read a fortran-formatted binary file that was saved at a defined precision (single or double).
#' @return  The read numeric matrix whose dimensions and precision are specified in the file
#' @param filename Name of the binary file to read
#' @param indexRow Vector of integers indicating the rows to be read from the file. Default \code{indexRow=NULL} will read all the rows
#' @param indexCol Vector of integers indicating the columns to be read from the file. Default \code{indexCol=NULL} will read all the columns
#' @param verbose \code{TRUE} or \code{FALSE} to whether printing file information
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

  # Read lines
  X <- .Call("readBinFile",filename,nsetRow,nsetCol,
      as.integer(indexRow),as.integer(indexCol))
  n <- X[[1]]; p <- X[[2]]; sizevar <- X[[3]];
  X <- X[[4]]
  if(verbose) cat("Loaded file '",filename,"'\n nrows=",n," ncols=",p," size=",sizevar," bytes\n")

  return(X)
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
  |  Sparse Family and Selection Index. Version 1.0.1                  |
  |====================================================================|
  ")
}
