
#====================================================================
#====================================================================
coef.SFI <- function(object,...)
{
  args0 <- list(...)

  df <- NULL
  if("df" %in% names(args0)) df <- args0$df

  if(!is.null(df)){
      if( df < 0 | df > length(object$trn) )
        stop("Parameter 'df' must be greater than zero and no greater than the number of elements in the training set")
  }
  if(ncol(object$df) == 1) df <- mean(as.vector(object$df))  # If only one lambda was considered

  which.df <- apply(object$df,1,function(x)which.min(abs(x-df)))
  if(is.character(object$BETA))
  {
    if(length(object$BETA)>1){
      index <- c(0,cumsum(object$nset))
    }else index <- c(0,length(object$tst))

    BETA <- c()
    for(i in seq_along(object$BETA))
    {
      if(!is.null(df))
      {
        indexRow <- which.df[seq(index[i]+1,index[i+1])]
        indexRow <- indexRow + c(0,cumsum(rep(ncol(object$df),length(object$tst)-1)))[1:length(indexRow)]
      }else indexRow <- NULL

      BETA <- rbind(BETA,readBinary(object$BETA[i],indexRow=indexRow,verbose=FALSE))
    }

    if(is.null(df)){
      BETA <- data.frame(BETA); colnames(BETA) <- NULL; rownames(BETA) <- NULL
      BETA <- split(BETA,rep(seq_along(object$tst),each=ncol(object$df)))
      BETA <- lapply(BETA,function(x) Matrix::Matrix(as.matrix(x), sparse=TRUE))

    }else BETA <- Matrix::Matrix(BETA, sparse=TRUE)

  }else{
      if(!is.null(df))
      {
        BETA <- matrix(NA,nrow=length(object$tst),ncol=length(object$trn))
        for(i in seq_along(object$BETA))  BETA[i,] <- object$BETA[[i]][which.df[i],]
        BETA <- Matrix::Matrix(BETA, sparse=TRUE)

      }else BETA <- object$BETA
  }

  BETA
}

#====================================================================
#====================================================================
fitted.SSI <- function(object,...)
{
  args0 <- list(...)
  if(length(args0)==0) stop("A matrix of predictors must be provided")

  X <- args0[[1]]
  yHat <- tcrossprod(X,as.matrix(object$beta))
  colnames(yHat) <- paste0("yHat",1:ncol(yHat))
  yHat
}

#====================================================================
#====================================================================
fitted.SFI <- function(object,...)
{
  args0 <- list(...)
  yTRN <- object$y[object$trn]-object$Xb[object$trn]
  index <- is.na(yTRN)
  if(length(index)>0) yTRN[index] <- 0

  if(is.character(object$BETA))
  {
    uHat <- c()
    for(i in seq_along(object$BETA))
    {
      uHat <- rbind(uHat,matrix(readBinary(object$BETA[i],verbose=FALSE) %*% yTRN,ncol=ncol(object$df),byrow=TRUE))
    }
  }else{
      uHat <- do.call("rbind",lapply(object$BETA,function(x)drop(as.matrix(x)%*%yTRN)))
  }
  dimnames(uHat) <- list(object$tst,paste0("yHat",1:ncol(uHat)))

  return(uHat)
}

#====================================================================
#====================================================================
plot.SFI <- function(...,title=NULL,py=c("accuracy","MSE"))
{
    PC1 <- PC2 <- PC1_TST <- PC2_TST <- PC1_TRN <- PC2_TRN <- loglambda <- NULL
    k <- model <- y <- trn_tst <- NULL
    py <- match.arg(py)
    args0 <- list(...)

    object <- args0[unlist(lapply(args0,function(x)class(x)=="SFI"))]
    if(length(object)==0) stop("No object of the class 'SFI' was provided")

    # Treat repeated fm$name
    modelNames <- unlist(lapply(object,function(x)x$name))
    index <- table(modelNames)[table(modelNames)>1]
    if(length(index)>0){
      for(i in seq_along(index)){
        tmp <- which(modelNames == names(index[i]))
        for(j in seq_along(tmp)) object[[tmp[j]]]$name <- paste0(object[[tmp[j]]]$name,".",j)
      }
    }

    trn <- do.call(rbind,lapply(object,function(x)x$trn))
    tst <- do.call(rbind,lapply(object,function(x)x$tst))
    if(any(apply(trn,2,function(x)length(unique(x)))!=1)) stop("'Training' set is not same across all 'SFI' objects")
    if(any(apply(tst,2,function(x)length(unique(x)))!=1)) stop("'Testing' set is not same across all 'SFI' objects")
    trn <- as.vector(trn[1,])
    tst <- as.vector(tst[1,])

    nTST <- length(tst)
    nTRN <- length(trn)

    theme0 <- theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.background = element_rect(fill="gray95"),
      #legend.key = element_rect(fill="gray95"),
      legend.box.spacing = unit(0.4, "lines"),
      legend.justification = c(1,ifelse(py=="MSE",1,0)),
      legend.position=c(0.96,ifelse(py=="MSE",0.96,0.04)),
      legend.title=element_blank(),
      legend.margin=margin(t=0,b=0.25,l=0.25,r=0.25,unit='line')
    )

      dat <- c()
      meanopt <- c()
      for(j in 1:length(object))
      {
          fm0 <- object[[j]]
          tmp <- summary.SFI(fm0)
          names(tmp[[py]]) <- paste0("SFI",1:length(tmp[[py]]))

          tt1 <- data.frame(SFI=names(tmp[[py]]),y=tmp[[py]],df=apply(fm0[['df']],2,mean),lambda=apply(fm0[['lambda']],2,mean))
          if(any(tt1$lambda<.Machine$double.eps)){
                tmp <- tt1$lambda[tt1$lambda>=.Machine$double.eps]
                tt1[tt1$lambda<.Machine$double.eps,'lambda'] <- ifelse(length(tmp)>0,min(tmp)/2,1E-6)
          }
          tt1 <- data.frame(obj=j,model=fm0$name,method=fm0$method,tt1,loglambda=-log(tt1$lambda),stringsAsFactors=FALSE)
          dat <- rbind(dat,tt1)

          index <- ifelse(py=="MSE",which.min(tt1$y),which.max(tt1$y))
          meanopt <- rbind(meanopt,tt1[index,])
      }

      dat$model <- factor(as.character(dat$model))

      # Models with a single SFI (G-BLUP or SFI with a given value of lambda)
      index <- unlist(lapply(split(dat,dat$obj),function(x)length(unique(x$SFI))))
      dat2 <- dat[dat$obj %in% names(index[index==1]),]

      # Remove data from models with a single SFI point
      dat <- dat[!dat$obj %in% names(index[index==1]),]
      meanopt <- meanopt[!meanopt$obj %in% names(index[index==1]),]
      dat <- dat[!(is.na(dat$y) | is.na(dat$loglambda)),]   # Remove NA values
      labY <- ifelse(py=="accuracy",expression('cor(y,'*hat(y)*')'),py)
      labX <- expression("-log("*lambda*")")

      if(nrow(dat)==0 | nrow(meanopt)==0)  stop("The plot can not be generated with the provided data")

      title0 <- paste0("Testing set ",py)
      if(!is.null(title)) title0 <- title

      # Labels and breaks for the DF axis
      ax2 <- getSecondAxis(dat$lambda,dat$df)
      brks0 <- ax2$breaks
      labs0 <- ax2$labels

      pt <- ggplot(dat,aes(loglambda,y,group=model,color=model)) + geom_line(size=0.66) +
          geom_hline(data=dat2,aes(yintercept=y,color=model,group=model),size=0.66) +
          labs(title=title0,x=labX,y=labY)+theme_bw() + theme0 + ylim(min(dat$y[dat$df>1]),max(dat$y)) +
          geom_vline(data=meanopt,aes(xintercept=loglambda),size=0.5,linetype="dotted",color="gray50")

      if(length(brks0)>3){
        pt <- pt + scale_x_continuous(sec.axis=sec_axis(~.+0,"Number of predictors",breaks=brks0,labels=labs0))
      }
      pt
}

#====================================================================
#====================================================================
plot.SFI_CV <- function(...,py=c("accuracy","MSE"), title=NULL,showFolds=FALSE)
{
    py <- match.arg(py)
    args0 <- list(...)
    obj_CV_fold <- negLogLambda <- CV <- fold <- y <- model <- NULL

    object <- args0[unlist(lapply(args0,function(x)class(x)=="SFI_CV"))]

    # Treat repeated fm$name
    modelNames <- unlist(lapply(object,function(x)unique(unlist(lapply(x,function(z)z$name)))))
    index <- table(modelNames)[table(modelNames)>1]
    if(length(index)>0){
      for(i in seq_along(index)){
        tmp <- which(modelNames == names(index[i]))
        for(j in seq_along(tmp)){
          for(k in 1:length(object[[tmp[j]]]))
          object[[tmp[j]]][[k]]$name <- paste0(object[[tmp[j]]][[k]]$name,".",j)
        }
      }
    }
    varNames <- c("y","df","lambda","negLogLambda")
    dat <- c()
    for(j in 1:length(object))
    {
        fm0 <- object[[j]]
        names(fm0) <- paste0("CV",1:length(fm0))
          rawdat <- do.call(rbind,lapply(fm0,function(x)reshape::melt(x[[py]])))
          colnames(rawdat) <- c("fold","SFI","y")
          rawdat <- data.frame(CV=unlist(lapply(strsplit(rownames(rawdat),"\\."),function(x)x[1])),rawdat)
          rawdat$df <- do.call(rbind,lapply(fm0,function(x)reshape::melt(x[['df']])))$value
          rawdat$lambda <- do.call(rbind,lapply(fm0,function(x)reshape::melt(x[['lambda']])))$value
          rawdat$negLogLambda <- -log(rawdat$lambda)
          rawdat <- data.frame(obj=j,model=fm0[[1]]$name,method=fm0[[1]]$method,rawdat,stringsAsFactors=FALSE)
        dat <- rbind(dat,rawdat)
    }

    # Average across folds
    avgdat <- do.call(rbind,lapply(split(dat,paste0(dat$obj,"_",dat$CV,"_",dat$SFI)),function(x){
        x[1,varNames] <- apply(as.matrix(x[,varNames]),2,mean,na.rm=TRUE)
        x$fold <- "mean"
        x[1,]
    }))

    # Average across partitions across folds
    overalldat <- do.call(rbind,lapply(split(avgdat,paste0(avgdat$obj,"_",avgdat$SFI)),function(x){
      x[1,varNames] <- apply(as.matrix(x[,varNames]),2,mean,na.rm=TRUE)
      x$CV <- "mean"
      x[1,]
    }))

    if(showFolds & length(object)==1){
      dat <- rbind(dat,avgdat,overalldat)
    }else dat <- rbind(avgdat,overalldat)

    dat$obj_CV_fold <- factor(paste0(dat$obj,"_",dat$CV,"_",dat$fold))
    dat$obj <- factor(as.character(dat$obj))

    # Optimum INDEX
    optdat <- do.call(rbind,lapply(split(overalldat,overalldat$obj),function(x){
      x[ifelse(py =="accuracy",which.max(x$y),which.min(x$y)),]
    }))
    optdat$obj_CV_fold <- factor(paste0(optdat$obj,"_",optdat$CV,"_",optdat$fold))

    # Models with a single SFI
    index <- unlist(lapply(split(dat,paste0(dat$obj,"_",dat$CV)),function(x)mean(table(x$fold))))
    dat2 <- dat[paste0(dat$obj,"_",dat$CV) %in% names(index[index == 1]),]

    # Remove data from models with a single SFI point
    dat <- dat[paste0(dat$obj,"_",dat$CV) %in% names(index[index > 1]),]
    optdat <- optdat[paste0(optdat$obj,"_",optdat$CV) %in% names(index[index > 1]),]

    dat <- dat[!is.na(dat$y) & !is.na(dat$negLogLambda),]   # Remove NA values
    labY <- ifelse(py=="accuracy",expression('cor(y,'*hat(y)*')'),py)
    labX <- expression("-log("*lambda*")")

    # Labels and breaks for the DF axis
    if(nrow(dat) == 0) stop("The plot cannot be generated with the provided data")
    ax2 <- getSecondAxis(dat$lambda,dat$df)
    brks0 <- ax2$breaks
    labs0 <- ax2$labels

    title0 <- paste0("Average cross-validated ",py)
    if(!is.null(title)) title0 <- title

    theme0 <- theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.background = element_rect(fill="gray95"),
        legend.box.spacing = unit(0.4, "lines"),
        legend.justification = c(1,ifelse(py=="MSE",1,0)),
        legend.position=c(0.97,ifelse(py=="MSE",0.97,0.03)),
        legend.title=element_blank(),
        legend.margin=margin(t=0,b=0.25,l=0.25,r=0.25,unit='line')
    )

    dat <- dat[dat$df > 0.99,]
    index1 <- which(dat$CV == "mean" & dat$fold == "mean")
    index2 <- which(dat$CV != "mean" & dat$fold == "mean")
    index3 <- which(dat$CV != "mean" & dat$fold != "mean")

    pt <- ggplot(dat,aes(group=obj_CV_fold)) +
             geom_hline(data=dat2[dat2$CV == "mean" & dat2$fold == "mean",],
                        aes(yintercept=y,color=model,group=obj_CV_fold),size=0.7) +
             labs(title=title0,x=labX,y=labY) + theme_bw() + theme0 +
             geom_vline(data=optdat,aes(xintercept=negLogLambda),size=0.5,linetype="dotted",color="gray50")

    if(showFolds){
      if(length(object)==1){
        pt <- pt + geom_line(data=dat[index3,],aes(negLogLambda,y,group=obj_CV_fold),color="gray70",size=0.3)
      }else cat("Results for individuals folds are not shown when plotting more than one model\n")
    }

    if(length(object)==1){   # Results from each CV
      pt <- pt + geom_line(data=dat[index2,],aes(negLogLambda,y,group=obj_CV_fold,color=CV),size=0.4) +
                 geom_line(data=dat[index1,],aes(negLogLambda,y,group=obj_CV_fold),color="gray5",size=0.7) +
                 theme(legend.position = "none")
    }else{
      pt <- pt + geom_line(data=dat[index1,],aes(negLogLambda,y,group=obj_CV_fold,color=model),size=0.7)
    }

    if(length(brks0)>3){
      pt <- pt + scale_x_continuous(sec.axis=sec_axis(~.+0,"Number of predictors",breaks=brks0,labels=labs0))
    }
    pt
}

#====================================================================
#====================================================================
summary.SFI_CV <- function(object, ...)
{
    args0 <- list(...)
    if(!inherits(object, "SFI_CV")) stop("The provided object is not from the class 'SFI'")

    df <- do.call(rbind,lapply(object,function(x)apply(x$df,2,mean,na.rm=TRUE)))
    lambda <- do.call(rbind,lapply(object,function(x)apply(x$lambda,2,mean,na.rm=TRUE)))
    MSE <- do.call(rbind,lapply(object,function(x)apply(x$MSE,2,mean,na.rm=TRUE)))
    accuracy <- do.call(rbind,lapply(object,function(x)apply(x$accuracy,2,mean,na.rm=TRUE)))
    rownames(df) <- rownames(lambda) <- rownames(accuracy) <- rownames(MSE) <- paste0("CV",1:nrow(df))

    out <- list(df=df,lambda=lambda,accuracy=accuracy,MSE=MSE)

    # Detect maximum accuracy by partition (curve)
    index <- apply(accuracy,1,which.max)
    tmp <- lapply(out,function(x)unlist(lapply(1:nrow(x),function(z)x[z,index[z]])))
    optCOR <- do.call(cbind,tmp)

    ##  Maximum accuracy averaging curves
    index <- which.max(apply(accuracy,2,mean,na.rm=TRUE))
    tmp <- unlist(lapply(out,function(x)apply(x,2,mean,na.rm=TRUE)[index]))
    optCOR <- rbind(optCOR,mean=tmp)

    # Detect minimum MSE by partition (curve)
    index <- apply(MSE,1,which.min)
    tmp <- lapply(out,function(x)unlist(lapply(1:nrow(x),function(z)x[z,index[z]])))
    optMSE <- do.call(cbind,tmp)

    ##  Minimum MSE averaging curves
    index <- which.min(apply(MSE,2,mean,na.rm=TRUE))
    tmp <- unlist(lapply(out,function(x)apply(x,2,mean,na.rm=TRUE)[index]))
    optMSE <- rbind(optMSE,mean=tmp)

    do.call(c, list(as.list(out), optCOR=list(optCOR), optMSE=list(optMSE)))
}

#====================================================================
#====================================================================
summary.SFI <- function(object,...)
{
    args0 <- list(...)
    if(!inherits(object, "SFI")) stop("The provided object is not from the class 'SFI'")

    tst <- object$tst
    y <- object$y

    df <- apply(object$df,2,mean)
    lambda <- apply(object$lambda,2,mean)
    uHat <- as.matrix(fitted.SFI(object))
    accuracy <- suppressWarnings(drop(stats::cor(y[tst],uHat,use="pairwise.complete.obs")))
    MSE <- suppressWarnings(apply((y[tst]-uHat)^2,2,sum,na.rm=TRUE)/length(tst))

    out <- data.frame(accuracy=accuracy,MSE=MSE,df=df,lambda=lambda)

    # Detect maximum accuracy
    index <- which.max(out$accuracy)
    if(length(index)==0)
    {
      optCOR <- out[1,]
      #if(nrow(out)>1) optCOR[1,] <- NA
    }else optCOR <- out[index,]

    # Detect minimum MSE
    index <- which.min(out$MSE)
    if(length(index)==0)
    {
      optMSE <- out[1,]
      #if(nrow(out)>1) optMSE[1,] <- NA
    }else optMSE <- out[index,]

    do.call(c, list(as.list(out), optCOR=list(optCOR), optMSE=list(optMSE)))
}
