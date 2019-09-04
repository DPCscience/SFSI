#' Sparse Family Index
#'
#' Computes the entire Elastic-Net solution for the regression coefficients of a Family Index simultaneously for all
#' values of the penalization parameter via either the Coordinate Descent (CD) or Least Angle Regression (LARS) algorithms.
#'
#' The model is fitted for each individual in the testing set which are predicted using all available observations of the response variable
#' from training data \eqn{\textbf{y}=(y_1,...,y_{nTRN})}. The coefficients are estimated as function of the the (genetic) variance among 
#' training data (\eqn{G_{TRN}}) and the (genetic) covariance between training individuals and each testing individual (\eqn{G_{TRN,TST(i)}}),
#' all taken from the genetic relatedness matrix (G).
#'
#' The \eqn{nTRN} coefficients \eqn{\boldsymbol{\beta}_i=(\beta_{i,1},...,\beta_{i,nTRN})} for the \eqn{i^{th}} individual are obtained by optimizing the function
#' \deqn{-G_{TRN,TST(i)}' \boldsymbol{\beta}_i + \frac{1}{2} \boldsymbol{\beta}_i'(G_{TRN} + ((1-h^2)/h^2)I) \boldsymbol{\beta}_i + \lambda J(\boldsymbol{\beta}_i)}
#' where \eqn{\lambda} is the penalization parameter and \eqn{J(\boldsymbol{\beta)}} is a penalty function given by
#' \deqn{\frac{1}{2}(1-\alpha)||\boldsymbol{\beta}_i||_2^2 + \alpha||\boldsymbol{\beta}_i||_1}
#' for \eqn{0\leq\alpha\leq 1}. Each individual solution is found using the 'SSI' function (see \code{help(SSI)} or \code{?SSI})
#' @return  List object containing the elements:
#' \itemize{
#'   \item \code{BETA}: list object containing, for each individual in testing set, a matrix of regression coefficients.
#'   \item \code{alpha}: value for the elastic-net weights used.
#'   \item \code{lambda}: matrix with the sequence of values of lambda used (for each individual in rows).
#'   \item \code{df}: degrees of freedom (averaged across individuals), number of non-zero predictors at each solution.
#'   \item \code{kernel}: transformation applied to the elements of \code{G}.
#' }
#' Elements used as inputs: \code{y}, \code{h2}, \code{training}, \code{testing}, \code{method}, \code{name}, are also returned. The returned object is of the class 'SFI' for which methods
#' \code{fitted}, \code{predict}, \code{plot} and \code{summary} exist
#' @param G Genetic relatedness matrix. This can be a name of a binary file where the matrix is storaged 
#' @param y Response variable
#' @param h2 Heritability of the response variable. Default is \code{h2=0.5}
#' @param training Vector of integers indicating which individuals are in training set. Default is \code{training=1:length(y)} will consider all individuals as training
#' @param testing Vector of integers indicating which individuals are in testing set. Default is \code{testing=1:length(y)} will consider all individuals as testing
#' @param indexG Vector of integers indicating which columns and rows will be read when \code{G} is the name of a binary file.
#' Default \code{indexG=NULL} will read the whole matrix
#' @param subset A two-elements numeric vector \eqn{c(j,J)} to fit the model only for the \eqn{j^{th}} subset out of \eqn{J} subsets that the
#' testing set will be divided into. Results can be automatically saved when \code{saveAt} parameter is provided and can be retrieved later
#' using function \code{collect}. Default is \code{subset=NULL} for no subsetting, then the model is fitted using all information
#' @param kernel Kernel transformation to be applied to \code{G[training,training]} . List consisting on one of:
#' \itemize{
#'   \item \code{list(kernel='GAUSSIAN',h)}. If \code{h} is not provided the value of \code{h=-2*log(0.5)} is used.
#'   \item \code{list(kernel='LAPLACIAN',h)}. If \code{h} is not provided the value of \code{h=-2*log(0.5)} is used.
#'   \item \code{list(kernel='POLYNOMIAL',a,b)}. The values of \code{a=1} and \code{b=2} are used when they are not provided.
#' }
#' Default \code{kernel=NULL} (no kernel)
#' @param method One of:
#' \itemize{
#'  \item \code{'CD1'}: Coordinate Descent algorithm that computes the coefficients for a provided grid of lambdas common to all individuals in testing set.
#'  \item \code{'CD2'}: Similar to \code{'CD1'} but using a grid of lambdas specific to each individual in testing set.
#'  \item \code{'LAR'}: Least Angle Regression algorithm that computes the entire sequence of all coefficients. Values of lambdas are calculated at each step.
#'  \item \code{'LAR-LASSO'}: Similar to \code{'LAR'} but solutions when a predictor leaves the solution are also returned.
#'  \item \code{'GBLUP'}: Coefficients are derived with no penalization and they correspond to those of the kinship-based BLUP
#' }
#' Default is \code{method='CD1'}
#' @param maxDF Maximum (average across individuals) number of predictors in the last solution (when \code{method='LAR'} or \code{'LAR-LASSO'}).
#' Default \code{maxDF=NULL} will calculate solutions including 1,2,...,nTRN predictors
#' @param lambda Penalization parameter sequence vector used for the Coordinate Descent algorithm.
#' Default is \code{lambda=NULL}, in this case a decreasing grid of \code{n='nLambda'} lambdas will be generated 
#' starting from a maximum equal to
#' \tabular{c}{\code{max(abs(G[training,testing])/alpha)}}
#' to a minimum equal to zero. If \code{alpha=0} the grid is generated starting from a maximum equal to 5. Only needed when \code{method='CD1'} or \code{'CD2'}
#' @param nLambda Number of lambdas generated when \code{lambda=NULL}
#' @param alpha Numeric between 0 and 1 indicating the weights for LASSO (alpha) and Ridge-Regression (1-alpha)
#' @param mc.cores Number of cores used to run the analysis in parallel. Default is \code{mc.cores=2}
#' @param tol Maximum error between two consecutive solutions of the iterative algorithm to declare convergence
#' @param maxIter Maximum number of iterations to run at each lambda step before convergence is reached
#' @param name Name given to the output for tagging purposes. Default \code{name=NULL} will give the name of the method used
#' @param saveAt Prefix name that will be added to the output files name to be saved, this may include a path. Regression coefficients
#' are saved as binary file as a single-precision (32 bits, 7 significant digits) variable. Default \code{saveAt=NULL} will no save any output
#' @param verbose \code{TRUE} or \code{FALSE} to whether printing each step
#' @examples
#' require(SFSI)
#' # Read data from BGLR package
#' data(wheat,package="BGLR")
#' X = scale(wheat.X)
#' G = tcrossprod(X)/ncol(X)    # Genomic relationship matrix
#' y = wheat.Y[,1]              # Response variable
#'
#' # Training and testing sets
#' set.seed(1234)
#' n = length(y)
#' pTST = 0.3      # percentage to predict
#' tst = sample(1:n,floor(pTST*n))
#' trn = (1:n)[-tst]
#'
#' # Calculate heritability
#' fm = BGLR::BGLR(y,ETA=list(list(K=G,model="RKHS")),nIter=10000,burnIn=3000,verbose=FALSE)
#' varU = fm$ETA[[1]]$varU
#' varE = fm$varE
#' h2 = varU/(varU + varE)
#'
#' # Sparse family index
#' fm = SFI(G,y,h2,trn,tst)
#'
#' yHat = fitted(fm)      # Predicted values or
#' pred = predict(fm)
#' yHat = pred$yHat       # Predicted values
#' corTST = cor(y[tst],yHat)[1,]  # Testing set accuracy
#' summary(fm)                    # Optimal accuracy
#' dat = data.frame(df=pred$df,logLambda=-log(pred$lambda),corTST)
#' plot(corTST ~ logLambda, data=dat[dat$df>1,]) # Accuracy along the regularization parameter
#' plot(fm)               # Same plot as above
#' plot(fm,py="MSE")      # MSE vs regularization
#' plot(fm,G=G)           # Coefficients path plot
#' plot(fm,G=G,PC=TRUE)   # Individual's training subsets
#' plot(fm,G=G,PC=TRUE,df=10)
#' @export
#' @references
#' \itemize{
#' \item \insertRef{Efron2004}{SFSI}
#' \item \insertRef{Friedman2007}{SFSI}
#' \item \insertRef{Hoerl1970}{SFSI}
#' \item \insertRef{Lush1947}{SFSI}
#' \item \insertRef{Perez2014}{SFSI}
#' \item \insertRef{Tibshirani1996}{SFSI}
#' \item \insertRef{VanRaden2008}{SFSI}
#' \item \insertRef{Zou2005}{SFSI}
#' }
#' @author Marco Lopez-Cruz (\email{lopezcru@msu.edu}) and Gustavo de los Campos
#' @importFrom parallel mclapply
#' @importFrom Rdpack reprompt
#' @importFrom Matrix Matrix
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords SFI

SFI <- function(G,y,h2=0.5,training=1:length(y),testing=1:length(y),indexG=NULL,subset=NULL,kernel=NULL,
    maxDF=NULL,lambda=NULL,nLambda=100,method=c("CD1","CD2","LAR","LAR-LASSO"),alpha=1,name=NULL,
    mc.cores=getOption("mc.cores", 2L),tol=1E-5,maxIter=800,saveAt=NULL,verbose=TRUE)
{
  method <- match.arg(method)

  nTRN <- length(training);  nTST <- length(testing)
  if(is.character(G)){
      G <- readBinary(G,indexRow=indexG,indexCol=indexG)
  }

  if(!is.matrix(G) | (length(G) != length(y)^2))
    stop("Object 'G' must be a squared matrix with number of rows (and columns) equal to the number of elements in 'y'")

  if(!is.double(G)) G <- apply(G,2,as.double)

  RHS <- G[training,testing,drop=FALSE]
  G <- G[training,training]

  for(i in 1:nTRN)  G[i,i] <- G[i,i] + (1-h2)/h2

  if(!is.null(kernel)){
    if(is.list(kernel) & is.null(kernel$kernel)) stop("Parameter 'kernel' must be a 'list' type object")
    G <- kernel2(G,kernel)
    kernel <- G$kernel
    G <- G$K
  }
  
  # Standardizing
  G <- scale_crossprod(G)
  sdx <- G$sdX
  G <- G$XtX
  RHS <- apply(RHS,2,function(x)x/sdx)

  if(is.null(lambda)){
    if(method == "CD1"){
        Cmax <- ifelse(alpha>0,max(abs(RHS)/alpha),5)
        lambda0 <- exp(seq(log(Cmax),log(1E-5),length=nLambda))
        lambda0[nLambda] <- 0
    }else lambda0 <- NULL
  }else{
    if(is.matrix(lambda)){
        if(nrow(lambda) != nTST) stop("Object 'lambda' must be a vector or a matrix with nTST rows")
    }else{
      lambda <- sapply(lambda,rep,times=nTST)
    }
  }

  if(is.null(name)) name <- paste(c(method,substr(kernel$kernel,1,2)),collapse="_")

  # Split the testing set into subsets. Only the subset provided will be fitted
  if(!is.null(subset)){
     if(!is.numeric(subset) & length(subset) != 2)
      stop("Object 'subset' must contain at least a 2-elements vector")
     partitions <- sort(rep(1:subset[2],ceiling(nTST/subset[2]))[1:nTST])
     index <- which(partitions == subset[1])
     tmp <- paste0(" of ",length(testing))
     testing <- testing[index]
     RHS <- RHS[,index,drop=FALSE]
  }else tmp <- ""

  cat(" Fitting SFI model for nTST=",length(testing),tmp," and nTRN=",length(training)," individuals\n",sep="")
  method <- ifelse(method %in% c("CD1","CD2"),"CD",method)

  compApply <- function(chunk)
  {
    rhs <- drop(RHS[,chunk])
    if(!is.null(lambda)) lambda0 <- lambda[chunk,]

    fm <- SSI(G,rhs,method=method,scale=FALSE,maxDF=maxDF,lambda=lambda0,
        nLambda=nLambda,alpha=alpha,tol=tol,maxIter=maxIter)

    # Returning betas to their original scale by scaling them by their SD
    B <- Matrix::Matrix(scale(fm$beta,FALSE,sdx), sparse=TRUE)

    cat(1,file=con,append=TRUE)
    if(verbose){
      utils::setTxtProgressBar(pb, nchar(scan(con,what="character",quiet=TRUE))/length(testing))
    }
    return(list(B=B,lambda=fm$lambda,testing=testing[chunk],df=fm$df))
  }

  pb = utils::txtProgressBar(style=3)
  con <- tempfile()
  if(mc.cores == 1L) {
    out = lapply(X=seq_along(testing),FUN=compApply)
  }else{
    out = parallel::mclapply(X=seq_along(testing),FUN=compApply,mc.cores=mc.cores)
  }
  close(pb); unlink(con)

  if(sum(testing != unlist(lapply(out,function(x)x$testing)))>0){
    stop("Matching error. Something went wrong during the analysis.")
  }

  index <- min(unlist(lapply(out,function(x)length(x$lambda))))
  out <- list(method=method, kernel=kernel, name=name, y=y, h2=h2,
              training=training,testing=testing, alpha=alpha,
              df=do.call("rbind",lapply(out,function(x)x$df[1:index])),
              lambda=do.call("rbind",lapply(out,function(x)x$lambda[1:index])),
              BETA=lapply(out,function(x) x$B[1:index, ,drop=FALSE])
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
