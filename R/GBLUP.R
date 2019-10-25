
GBLUP <- function(G,y,h2=0.5,trn=1:length(y),tst=1:length(y),kernel=NULL)
{
	yTRN <- y[trn]-mean(y[trn],na.rm=TRUE)
	index <- is.na(yTRN)
	if(length(index)>0) yTRN[index] <- 0

	## Coeff. matrix
	RHS <- G[trn,tst]
	G <- G[trn,trn]

	for(i in 1:length(trn))  G[i,i] <- G[i,i]+(1-h2)/h2

	if(!is.null(kernel)){
	  if(is.list(kernel) & is.null(kernel$kernel))  stop("Parameter 'kernel' must be a 'list' type object")
	  kernel2(G,kernel,void=TRUE)
	}

	Ginv <- chol2inv(chol(G))
  BETA <- crossprod(RHS,Ginv)

	yHat <- drop(BETA%*%yTRN)
	correlation <- cor(yHat,y[tst])
	MSE <- sum((yHat-y[tst])^2)/length(tst)

  out <- list(y=y, h2=h2, trn=trn, tst=tst, BETA=BETA,
    yHat=yHat, correlation=correlation, accuracy=correlation/sqrt(h2),
    MSE=MSE, kernel=kernel, method="GBLUP")
  class(out) <- "SFI"
	return(out)
}
