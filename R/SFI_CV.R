#' Sparse Family Index Cross-Validation
#'
#' Fit the Sparse Family Index in a cross validation fashion by randomly splitting data into 2,3,4,5, or 10 folds.
#' These partitions are equivalent to 50:50, 66:33, 75:25, 80:20, and 90:10 % training-testing
#' @return  List object containing the elements:
#' \itemize{
#'   \item folds: matrix containing the folds used for the cross-validation.
#'   \item correlation: matrix with the correlation between observed and predicted values (in testing set) within each fold (in rows).
#'   \item accuracy: matrix with the accuracy (in testing set) within each fold (in rows).
#'   \item MSE: matrix with the mean squared error of prediction (in testing set) within each fold (in rows).
#'   \item lambda: matrix with the sequence of values of lambda used (averaged across individuals) within each fold (in rows).
#'   \item df: matrix with the degrees of freedom (averaged across individuals) within each fold (in rows).
#'   \item kernel: transformation applied to the elements of 'G'.
#' }
#' Elements used as inputs: 'h2','method','name', are also returned. The returned object is of the class 'SFI_CV' for which methods 'plot' and 'summary' exist
#' @param G Genetic relatedness matrix
#' @param y Response variable
#' @param h2 Heritability of the response variable
#' @param training Index for the individuals for which the cross-validation will be implemented.
#' Default is \eqn{training=1:length(y)} will consider all individuals
#' @param nFolds Number of non-overlaping folds in which the data is splitted. Options 2,3,5,10
#' @param kernel Kernel transformation applied to 'G'. List consisting on one of:
#' \itemize{
#'   \item list(kernel='GAUSSIAN',h). If \eqn{h} is not provided the value of \eqn{h=-2*log(0.5)} is used.
#'   \item list(kernel='LAPLACIAN',h). If \eqn{h} is not provided the value of \eqn{h=-2*log(0.5)} is used.
#'   \item list(kernel='POLYNOMIAL',a,b). The values of \eqn{a=1} and \eqn{b=2} are used when they are not provided.
#' }
#' Default \eqn{kernel=NULL} (no kernel)
#' @param method One of:
#' \itemize{
#'  \item 'CD1': Coordinate Descent algorithm that computes the coefficients for a provided grid of lambdas common to all individuals in testing set.
#'  \item 'CD2': Similar to 'CD1' but using a grid of lambdas specific to each individual in testing set.
#'  \item 'LAR': Least Angle Regression algorithm that computes the entire sequence of all coefficients. Values of lambdas are calculated at each step.
#'  \item 'LAR-LASSO': Similar to 'LAR' but solutions when a predictor leaves the solution are also returned.
#'  \item 'GBLUP': Coefficients are derived with no penalization and they correspond to those of the kinship-based-BLUP
#' }
#' @param maxDF Maximum (average across individuals) number of predictors in the last solution (when method='LAR' or 'LAR-LASSO').
#' Default \eqn{maxDF=NULL} will calculate solutions including 1,2,...,nTRN predictors
#' @param lambda Penalization parameter sequence vector. Default is \eqn{lambda=NULL}, in this case a decreasing grid of
#' n='nLambda' lambdas will be generated starting from a maximum equal to \eqn{max(abs(G[training,testing])/alpha)} to a minumum equal to zero.
#' If \eqn{alpha=0} the grid is generated starting from a maximum equal to 5.
#' Only needed for method='CD1' or 'CD2'
#' @param nLambda Number of lambdas generated when \eqn{lambda=NULL}
#' @param alpha Numeric between 0 and 1 indicating the weights for LASSO (alpha) and Ridge-Regression (1-alpha)
#' @param nCores Number of cores used to run the analysis in parallel. Default is \eqn{nCores=2}
#' @param tol Maximum error between two consecutive solutions of the iterative algorithm to declare convergence
#' @param maxIter Maximum number of iterations to run at each lambda step before convergence is reached
#' @param seed Numeric seed to fix randomization when creating folds
#' @param name Name given to the output for tagging purposes. Default \eqn{name=NULL} will give the name of the method used
#' @param verbose TRUE or FALSE to whether printing each step
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
#' # Obtain a value of lambda via cross validation in training set
#' fm1 = SFI_CV(G,y,h2,trn)
#' plot(fm1)
#' lambda = summary(fm1)[[1]][[1]]$max[1,'lambda']
#' fm2 = SFI(G,y,h2,trn,tst,lambda=lambda)
#' predict(fm2)$correlation        # Testing set accuracy
#'
#' fm3 = SFI(G,y,h2,trn,tst,lambda=0)
#' predict(fm3)$correlation        # Testing set accuracy of the non-sparse index
#' @export
#' @references
#' \itemize{
#' \item \insertRef{Efron2004}{SFSI}
#' \item \insertRef{Friedman2007}{SFSI}
#' \item \insertRef{Lush1947}{SFSI}
#' \item \insertRef{Perez2014}{SFSI}
#' \item \insertRef{VanRaden2008}{SFSI}
#' }
#' @author Marco Lopez-Cruz (\email{lopezcru@@msu.edu}) and Gustavo de los Campos
#' @keywords SFI_CV

SFI_CV <- function(G,y,h2=0.5,training=1:length(y),nFolds=5,indexG=NULL,kernel=NULL,maxDF=NULL,
    lambda=NULL,nLambda=100,method=c("CD1","CD2","LAR","LAR-LASSO","GBLUP"),alpha=1,
    nCores=getOption("mc.cores", 2L),tol=2E-5,maxIter=750,seed=123,name=NULL,verbose=TRUE)
{
    nFolds <- match.arg(as.character(nFolds),choices=c(2,3,4,5,10))
    method <- match.arg(method)

    if(is.character(G)){
        G <- readBinary(G,indexRow=indexG,indexCol=indexG)
    }

    if(!is.null(kernel)){
        if(is.list(kernel) & is.null(kernel$kernel)) stop("Parameter 'kernel' must be a 'list' type object")
        G <- kernel2(G,kernel)
        kernel <- G$kernel
        G <- G$K
    }

    if(is.null(name)) name <- paste(c(method,substr(kernel$kernel,1,2)),collapse="_")
    nFolds <- as.numeric(nFolds)
    nTRN <- length(training)

    # Creating folds
    set.seed(seed)
    folds <- rep(seq(1:nFolds),ceiling(nTRN/nFolds))[1:nTRN]
    folds <- sample(folds)

    out <- vector("list",nFolds)
    for(j in 1:nFolds)
    {
        if(verbose) cat(" Fold ",j," of ",nFolds,"\n",sep="")
        if(method == "GBLUP"){
            fm <- GBLUP(G,y,h2,training[folds!=j],training[folds==j])
        }else{
            fm <- SFI(G,y,h2,training[folds!=j],training[folds==j],maxDF=maxDF,nCores=nCores,alpha=alpha,
                method=method,lambda=lambda,nLambda=nLambda,tol=tol,maxIter=maxIter,verbose=verbose)
        }
        fv <- predict(fm)
        out[[j]] <- list(correlation=fv$correlation, accuracy=fv$accuracy, MSE=fv$MSE,df=fv$df,lambda=fv$lambda)
    }
    index <- min(unlist(lapply(out,function(x)length(x$lambda))))
    out <- list(method=method, kernel=kernel, name=name,y=y, h2=h2,training=training,
                alpha=alpha, folds=data.frame(training,fold=folds),
                correlation=do.call("rbind",lapply(out,function(x)x$correlation[1:index])),
                accuracy=do.call("rbind",lapply(out,function(x)x$accuracy[1:index])),
                MSE=do.call("rbind",lapply(out,function(x)x$MSE[1:index])),
                df=do.call("rbind",lapply(out,function(x)x$df[1:index])),
                lambda=do.call("rbind",lapply(out,function(x)x$lambda[1:index]))
              )
    class(out) <- "SFI_CV"
    return(out)
}
