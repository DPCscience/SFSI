#' GBLUP function
#'
#' Fits the kinship-based Best Linear Unbiased Prediction (BLUP) model.
#'
#' @param G Genetic relatedness matrix
#' @param y Response variable
#' @param h2 Heritability of the response variable. Default is \code{h2=0.5}
#' @param training Vector of integers indicating which individuals are in training set. Default is \code{training=1:length(y)} will consider all individuals as training
#' @param testing Vector of integers indicating which individuals are in testing set. Default is \code{testing=1:length(y)} will consider all individuals as testing
#' @param kernel Kernel transformation to be applied to \code{G[training,training]}. List consisting on one of:
#' \itemize{
#'   \item \code{list(kernel='GAUSSIAN',h)}. If \code{h} is not provided the value of \code{h=-2*log(0.5)} is used.
#'   \item \code{list(kernel='LAPLACIAN',h)}. If \code{h} is not provided the value of \code{h=-2*log(0.5)} is used.
#'   \item \code{list(kernel='POLYNOMIAL',a,b)}. The values of \code{a=1} and \code{b=2} are used when they are not provided.
#' }
#' Default \code{kernel=NULL} (no kernel)
#' @examples
#' set.seed(1234)
#' require(SFSI)
#' # Read data from BGLR package
#' data(wheat,package="BGLR")
#' X = scale(wheat.X)
#' G = tcrossprod(X)/ncol(X)    # Genomic relationship matrix
#' y = wheat.Y[,1]              # Response variable
#'
#' # Training and testing sets
#' set.seed(1234)
#' tst = sample(seq_along(y),150)     # 150 individuals to predict
#' trn = seq_along(y)[-tst]           # 449 individuals to train model
#'
#' # Calculate heritability
#' fm = BGLR::BGLR(y,ETA=list(list(K=G,model="RKHS")),nIter=10000,burnIn=3000,verbose=FALSE)
#' varU = fm$ETA[[1]]$varU
#' varE = fm$varE
#' h2 = varU/(varU + varE)
#'
#' fm = GBLUP(G,y,h2,trn,tst)
#' plot(y[tst],fitted(fm))        # Predicted vs observed values in testing set
#' predict(fm)$correlation        # Testing set accuracy
#' @keywords BGLUP
#' @export
#' @references
#' \itemize{
#' \item \insertRef{Perez2014}{SFSI}
#' \item \insertRef{VanRaden2008}{SFSI}
#' }
#' @author Marco Lopez-Cruz (\email{lopezcru@msu.edu}) and Gustavo de los Campos

GBLUP <- function(G,y,h2,training=1:length(y),testing=1:length(y),kernel=NULL)
{
	yTRN <- y[training]-mean(y[training])

	## Coeff. matrix
	RHS <- G[training,testing]
	G <- G[training,training]

	for(i in 1:length(training))  G[i,i] <- G[i,i]+(1-h2)/h2
	
	if(!is.null(kernel)){
	  if(is.list(kernel) & is.null(kernel$kernel))  stop("Parameter 'kernel' must be a 'list' type object")
	  G <- kernel2(G,kernel)
	  kernel <- G$kernel
	  G <- G$K
	}

	Ginv <- chol2inv(chol(G))
  BETA <- crossprod(RHS,Ginv)

	yHat <- drop(BETA%*%yTRN)
	correlation <- cor(yHat,y[testing])
	MSE <- sum((yHat-y[testing])^2)/length(testing)

  out <- list(y=y, h2=h2, training=training, testing=testing, BETA=BETA,
    yHat=yHat, correlation=correlation, accuracy=correlation/sqrt(h2),
    MSE=MSE, kernel=kernel, method="GBLUP")
  class(out) <- "SFI"
	return(out)
}
