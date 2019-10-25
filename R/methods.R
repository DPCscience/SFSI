
#====================================================================
#====================================================================
coef.SFI <- function(object,...)
{
  args0 <- list(...)
  
  df <- NULL
  if("df" %in% names(args0)) df <- args0$df
  
  if(object$method=="GBLUP")
  {
      BETA <- Matrix::Matrix(object$BETA, sparse=TRUE)
  }else{
    if(!is.null(df)){
      if( df < 0 | df > length(object$trn) )
        stop("Parameter 'df' must be greater than zero and no greater than the number of elements in the training set")
    }
    df0 <- apply(object$df,1,function(x)which.min(abs(x-df)))
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
          indexRow <- df0[seq(index[i]+1,index[i+1])]
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
        BETA <- c()
        for(i in seq_along(object$BETA))  BETA <- rbind(BETA,object$BETA[[i]][df0[i],])
        BETA <- Matrix::Matrix(BETA, sparse=TRUE)

      }else BETA <- object$BETA
    }
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
  colnames(yHat) <- paste0("SSI",1:ncol(yHat))
  yHat
}

#====================================================================
#====================================================================
fitted.SFI <- function(object,...)
{
  args0 <- list(...)
  yTRN <- object$y[object$trn]-mean(object$y[object$trn],na.rm=TRUE)
  index <- is.na(yTRN)
  if(length(index)>0) yTRN[index] <- 0

  if(object$method == "GBLUP"){
    yHat <- drop(object$BETA%*%yTRN)
  }else{
    if(is.character(object$BETA)){
      yHat <- c()
      for(i in seq_along(object$BETA))
      {
        yHat <- rbind(yHat,matrix(readBinary(object$BETA[i],verbose=FALSE) %*% yTRN,ncol=ncol(object$df),byrow=TRUE))
      }
    }else{
       yHat <- do.call("rbind",lapply(object$BETA,function(x)drop(as.matrix(x)%*%yTRN)))
    }
    dimnames(yHat) <- list(object$tst,paste0("SFI",1:ncol(yHat)))
  }
  return(yHat)
}

#====================================================================
#====================================================================
plot.SFI <- function(...,df=NULL,G=NULL,path=FALSE,title=NULL,maxCor=0.85,
  py=c("correlation","accuracy","MSE"))
{
    PC1 <- PC2 <- PC1_TST <- PC2_TST <- PC1_TRN <- PC2_TRN <- NULL
    g <- model <- y <- trn_tst <- NULL
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

    flagPC <- ifelse(!is.null(G),!path,FALSE)
    flagPath <- path
    flagCor <- is.null(G) & !path

    if(flagPC | flagPath){
      if(length(object)>1) cat("Only the fist 'SFI' object is considered\n")
      object <- object[[1]]
      if(!is.null(G) & !is.null(object$kernel)) G <- kernel2(G,kernel=object$kernel)
    }

    if(flagPC){
        PC <-  RSpectra::eigs_sym(G, 2)
        expvarPC <- 100*(PC$values^2)/sum(G^2)
        labelsPC <- paste0("PC",1:2,paste0(" (",sprintf('%.1f',expvarPC[1:2]),"%)"))
    }

    theme0 <- theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.background = element_rect(fill="gray95"),
        #legend.key = element_rect(fill="gray95"),
        legend.box.spacing = unit(0.4, "lines")
    )

    if(flagCor){
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
      labY <- ifelse(py=="correlation",expression('cor(y,'*hat(y)*')'),py)
      labX <- expression("-log("*lambda*")")

      if(nrow(dat)==0 | nrow(meanopt)==0)  stop("The plot can not be generated with the provided data")

      title0 <- paste0("Testing set ",py)
      if(!is.null(title)) title0 <- title

      # Labels and breaks for the DF axis
      breaks0 <- seq(min(dat$loglambda),max(dat$loglambda),length=6)
      labels0 <- floor(stats::predict(stats::smooth.spline(dat$loglambda,dat$df),breaks0)$y)
      labels0[1] <- 1
      breaks0 <- stats::predict(stats::smooth.spline(dat$df, dat$loglambda),labels0)$y

      theme0 <- theme0 + theme(legend.justification = c(1,ifelse(py=="MSE",1,0)),
                  legend.position=c(0.96,ifelse(py=="MSE",0.96,0.04)),
                  legend.title=element_blank(),
                  legend.margin=margin(t=0,b=0.25,l=0.25,r=0.25,unit='line'))

      pt <- ggplot(dat[dat$df>1,],aes(loglambda,y,group=model,color=model)) + geom_line(size=0.66) +
          geom_hline(data=dat2,aes(yintercept=y,color=model,group=model),size=0.66) +
          labs(title=title0,x=labX,y=labY)+theme_bw() + theme0 +
          geom_vline(data=meanopt,aes(xintercept=loglambda),size=0.5,linetype="dotted",color="gray50")+
          scale_x_continuous(sec.axis = sec_axis(~.+0,"number of predictors",breaks=breaks0,labels=labels0))
      print(pt)
    }

    if(flagPC)    # Network plot
    {
        if(is.null(df)) df <- summary.SFI(object)$opt$df

        BETA <- coef.SFI(object,df=df)
        df <- mean(apply(BETA,1,function(y)sum(abs(y)>0)))

        PC <- PC$vectors
        dat <- data.frame(id=1:nrow(PC),PC[,1:2])

        set <- ifelse((1:nrow(dat)) %in% tst, 'testing','inactive')
        index <- which(apply(BETA,2,function(x)any(abs(x)>0)))   # Active set
        set[trn[index]] <- "active"
        dat <- data.frame(dat,set)
        colnames(dat) <- c("id","PC1","PC2","set")

        setLevels <- c("inactive","active","testing")
        dat$set <- factor(as.character(dat$set),levels=setLevels)
        index <- setLevels %in% unique(as.character(dat$set))

        title0 <- bquote(.(object$name)*". Number of predictors="*.(round(df)))
        if(!is.null(title)) title0 <- title

        pt <- ggplot(dat,aes(x=PC1,y=PC2)) +
          labs(title=title0,x=labelsPC[1],y=labelsPC[2])+
          geom_point(data=dat[dat$set=="inactive",],aes(fill=set),colour="gray35",shape=21,size=1)

        for(i in 1:nTST)
        {
            indexTRN <- which(abs(BETA[i,])>0)
            if(length(indexTRN)>0)
            {
                dat1 <- dat[trn,c("PC1","PC2")][indexTRN,]
                dat2 <- dat[tst,c("PC1","PC2")][i,]
                colnames(dat1) <- paste0(colnames(dat1),"_TRN")
                colnames(dat2) <- paste0(colnames(dat2),"_TST")
                dat1 <- data.frame(dat2[rep(1,nrow(dat1)),],dat1)
                dat1$effect <- abs(BETA[i,indexTRN])
                dat1$effect <-  dat1$effect /max(dat1$effect)
                pt <- pt + geom_segment(data=dat1,aes(x=PC1_TST,y=PC2_TST,xend=PC1_TRN,yend=PC2_TRN),
                    size=0.15,color="gray65") 
            }
        }

        pt <- pt  +
            geom_point(data=dat[dat$set=="active",],aes(fill=set),shape=21,colour="black",size=2) +
            geom_point(data=dat[dat$set=="testing",],aes(fill=set),shape=21,colour="black",size=3) +
            scale_fill_manual(values=c("inactive"="gray35","active"="#56B4E9","testing"="#E69F00")[index]) + theme_bw() +
            theme0 + theme(legend.title=element_blank(),legend.margin=margin(t=0,b=0.25,l=0.25,r=0.25,unit='line'))
        print(pt)

    }

    if(flagPath){  # Coefficients path plot
        BETA <- coef.SFI(object)
        lambda <- apply(object$lambda,2,mean)
        df <- apply(object$df,2,mean)
        nDF <- length(df)
        if(lambda[nDF] < .Machine$double.eps*1000)  lambda[nDF] <- lambda[nDF-1]/2
        if(nDF==1) stop("Coefficients path plot can not be generated for 'nLambda=1'")

        dat <- c()  # Beta coefficients path
        for(i in seq_along(BETA))
        {
            B0 <- as.matrix(BETA[[i]])
            index <- getIndexCorrelated(B0,maxCor)

            nTRN0 <- length(index)
            if(nTRN0 > 1){
                g0 <- NA
                if(!is.null(G)) g0 <- rep(G[tst[i],trn[index]],each=nDF)
                tst0 <- factor(rep(tst[i],nTRN0*nDF))
                trn0 <- factor(rep(trn[index],each=nDF))
                trn_tst0 <- tst0:trn0
                df0 <- rep(df,nTRN0)
                lambda0 <- rep(lambda,nTRN0)
                B0 <- as.vector(B0[,index])
                dat <- rbind(dat,data.frame(df=df0,lambda=lambda0,beta=B0,g=g0,tst=tst0,trn_tst=trn_tst0))
            }
        }

        # Labels and breaks for the DF axis
        loglambda <- -log(lambda)
        breaks0 <- seq(min(loglambda),max(loglambda),length=6)
        labels0 <- floor(stats::predict(stats::smooth.spline(loglambda,df),breaks0)$y)
        labels0[1] <- 1
        breaks0 <- stats::predict(stats::smooth.spline(df, loglambda),labels0)$y

        title0 <- bquote("Coefficients path. "*.(object$name))
        if(!is.null(title)) title0 <- title

        theme0 <- theme0 + theme(legend.key.height=unit(3,"line"),
                      legend.key.width = unit(0.8, "lines"))

        if(!is.null(G))
        {
          pt <- ggplot(dat,aes(-log(lambda),beta,color=g,group=trn_tst))+ viridis::scale_color_viridis() +
              geom_line() + theme_bw() + theme0 +
              labs(title=title0,y=expression(beta),x=expression("-log("*lambda*")")) +
              scale_x_continuous(sec.axis = sec_axis(~.+0,"number of predictors",
                breaks=breaks0,labels=labels0))
        }else{
          pt <- ggplot(dat,aes(-log(lambda),beta,color=trn_tst,group=trn_tst))+
            geom_line() + theme_bw() + theme0 +
            labs(title=title0,y=expression(beta),x=expression("-log("*lambda*")")) +
            scale_x_continuous(sec.axis = sec_axis(~.+0,"number of predictors",
                                                   breaks=breaks0,labels=labels0)) +
            theme(legend.position = "none")
        }
        print(pt)
    }
}

#====================================================================
#====================================================================
plot.SFI_CV <- function(...,py=c("correlation","accuracy","MSE"),
  title=NULL,showFolds=FALSE)
{
    py <- match.arg(py)
    args0 <- list(...)
    fold <- loglambda <- y <- model <- NULL

    object <- args0[unlist(lapply(args0,function(x)class(x)=="SFI_CV"))]

    nFolds <- unlist(lapply(object,function(x)length(unique(x$folds$fold))))
    if(length(unique(nFolds))>1) stop("The number of folds are different across all 'SFI_CV' objects")
    nFolds <- nFolds[1]
    colors <- c(scales::hue_pal()(nFolds),"black")

    # Treat repeated fm$name
    modelNames <- unlist(lapply(object,function(x)x$name))
    index <- table(modelNames)[table(modelNames)>1]
    if(length(index)>0){
      for(i in seq_along(index)){
        tmp <- which(modelNames == names(index[i]))
        for(j in seq_along(tmp)) object[[tmp[j]]]$name <- paste0(object[[tmp[j]]]$name,".",j)
      }
    }

    dat <- c()
    meanopt <- c()
    for(j in 1:length(object))
    {
        fm0 <- object[[j]]
        colnames(fm0[[py]]) <- paste0("SFI",1:ncol(fm0[[py]]))

        tt1 <- reshape::melt(fm0[[py]])
        tt1$X2 <- as.character(tt1$X2)
        meany <- apply(fm0[[py]],2,mean)
        if(fm0$method=="GBLUP"){
            tt1$df <- NA
            tt1$lambda <- NA
            meandf <- meanlambda <- NA
        }else{
            if(any(fm0[['lambda']]<.Machine$double.eps)){
                  tmp <- fm0[['lambda']][fm0[['lambda']]>=.Machine$double.eps]
                  fm0[['lambda']][fm0[['lambda']]<.Machine$double.eps] <- ifelse(length(tmp)>0,min(tmp)/2,1E-6)
            }
            tt1$df <- reshape::melt(fm0[['df']])$value
            tt1$lambda <- reshape::melt(fm0[['lambda']])$value
            meandf <- apply(fm0[['df']],2,mean)
            meanlambda <- apply(fm0[['lambda']],2,mean)
        }
        tt1 <- data.frame(obj=j,model=fm0$name,method=fm0$method,tt1,loglambda=-log(tt1$lambda),stringsAsFactors=FALSE)
        tt2 <- data.frame(obj=j,model=fm0$name,method=fm0$method,X1="mean",X2=names(meany),
               value=meany,df=meandf,lambda=meanlambda,loglambda=-log(meanlambda),stringsAsFactors=FALSE)
        dat <- rbind(dat,tt1,tt2)

        index <- ifelse(py=="MSE" & fm0$method!="GBLUP",which.min(tt2$value),which.max(tt2$value))
        meanopt <- rbind(meanopt,tt2[index,])
    }
    names(dat) <- names(meanopt) <-  c("obj","model","method","fold","SFI","y","df","lambda","loglambda")
    dat$fold <- factor(as.character(dat$fold),levels=c(1:nFolds,"mean"))

    # Order factor models to have G-BLUP first (if exists)
    tmp <- do.call("rbind",lapply(split(dat,dat$obj),function(x)x[1,]))
    tmp <- tmp[order(tmp$method %in% "GBLUP"),]
    dat$model <- factor(as.character(dat$model),levels=as.character(tmp$model))

    # Models with a single SFI (G-BLUP or SFI with a given value of lambda)
    index <- unlist(lapply(split(dat,dat$obj),function(x)mean(table(x$fold))))
    dat2 <- dat[dat$obj %in% names(index[index==1]),]

    # Remove data from models with a single SFI point
    dat <- dat[!dat$obj %in% names(index[index==1]),]
    meanopt <- meanopt[!meanopt$obj %in% names(index[index==1]),]
    dat <- dat[!(is.na(dat$y) | is.na(dat$loglambda)),]   # Remove NA values
    labY <- ifelse(py=="correlation",expression('cor(y,'*hat(y)*')'),py)
    labX <- expression("-log("*lambda*")")

    # Labels and breaks for the DF axis
    breaks0 <- seq(min(dat$loglambda),max(dat$loglambda),length=6)
    labels0 <- floor(stats::predict(stats::smooth.spline(dat$loglambda,dat$df),breaks0)$y)
    labels0[1] <- 1
    breaks0 <- stats::predict(stats::smooth.spline(dat$df, dat$loglambda),labels0)$y

    title0 <- paste0(nFolds," folds CV. ",ifelse(!showFolds,"Average testing","Testing")," set ",py)
    if(!is.null(title)) title0 <- title

    theme0 <- theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.background = element_rect(fill="gray95"),
        legend.box.spacing = unit(0.4, "lines")
    )
    dat <- dat[dat$df>1,]
    if(showFolds)
    {
        pt <- ggplot(dat,aes(loglambda,y,group=fold:model,color=fold,linetype=model)) + geom_line(size=0.66) +
            geom_hline(data=dat2,aes(yintercept=y,color=fold,group=fold:model,linetype=model),size=0.66) +
            labs(title=title0,x=labX,y=labY)+theme_bw() + scale_color_manual(breaks=levels(dat$fold),values=colors)+
            geom_vline(data=meanopt,aes(xintercept=loglambda),size=0.5,linetype="dotted",color="gray50")+theme0 +
            scale_x_continuous(sec.axis = sec_axis(~.+0,"number of predictors",breaks=breaks0,labels=labels0))

    }else{
        dat <- dat[dat$fold=="mean",]
        dat2 <- dat2[dat2$fold == "mean",]
        theme0 <- theme0 + theme(legend.justification = c(1,ifelse(py=="MSE",1,0)),
                    legend.position=c(0.96,ifelse(py=="MSE",0.96,0.04)),
                    legend.title=element_blank(),
                    legend.margin=margin(t=0,b=0.25,l=0.25,r=0.25,unit='line'))

        pt <- ggplot(dat,aes(loglambda,y,group=model,color=model)) + geom_line(size=0.66) +
            geom_hline(data=dat2,aes(yintercept=y,group=model,color=model),size=0.66) +
            labs(title=title0,x=labX,y=labY)+theme_bw() +
            scale_color_manual(breaks=levels(dat$model),values=scales::hue_pal()(nlevels(dat$model)))+
            geom_vline(data=meanopt,aes(xintercept=loglambda,color=model),size=0.5,linetype="dotted")+theme0+
            scale_x_continuous(sec.axis = sec_axis(~.+0,"number of predictors",breaks=breaks0,labels=labels0))
    }
    print(pt)
}

#====================================================================
#====================================================================
plot.SSI <- function(x,...)
{
    if(!inherits(x, "SSI")) stop("The provided object is not from the class 'SSI'")

    theme0 <- theme(
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none"
    )
    args0 <- list(...)

    beta <- as.vector(x$beta)
    id <- factor(rep(seq(ncol(x$beta)),each=nrow(x$beta)))
    df <- x$df
    lambda <- x$lambda
    if(lambda[length(lambda)] < .Machine$double.eps*1000) lambda[length(lambda)] <- lambda[length(lambda)-1]/2

    dat <- data.frame(df=rep(df,ncol(x$beta)),lambda=rep(lambda,ncol(x$beta)),beta=beta,id=id)

    # Labels and breaks for the DF axis
    loglambda <- -log(lambda)
    breaks0 <- seq(min(loglambda),max(loglambda),length=6)
    labels0 <- floor(stats::predict(stats::smooth.spline(loglambda,df),breaks0)$y)
    labels0[1] <- 1
    breaks0 <- stats::predict(stats::smooth.spline(df, loglambda),labels0)$y

    title0 <- bquote("Coefficients path. "*.(x$name))
    if("title" %in% names(args0)) title0 <- args0$title
    
    pt <- ggplot(dat,aes(-log(lambda),beta)) +
        scale_x_continuous(sec.axis = sec_axis(~.+0,"number of predictors",
                breaks=breaks0,labels=labels0))+
        geom_line(size=0.5,aes(color=id,group=id)) + theme_bw() + theme0 +
        labs(title=title0,y=expression(beta),x=expression("-log("*lambda*")"))
    print(pt)
}

#====================================================================
#====================================================================
summary.SFI_CV <- function(object, ...)
{
    args0 <- list(...)
    if(!inherits(object, "SFI_CV")) stop("The provided object is not from the class 'SFI'")
    
    if(object$method == "GBLUP"){
      df <- NA; lambda <- NA
    }else{
      df <- apply(object$df,2,mean)
      lambda <- apply(object$lambda,2,mean)
    }
    
    correlation <- apply(object$correlation,2,mean)
    MSE <- apply(object$MSE,2,mean)
    accuracy <- apply(object$accuracy,2,mean)
    
    out <- data.frame(correlation=correlation,accuracy=accuracy,MSE=MSE,df=df,lambda=lambda)
    
    # Detect maximum correlation
    index <- which.max(out$correlation)
    if(length(index)==0)
    {
      if(nrow(out)==1){
        index <- 1
      }else stop("'summary' method can not be implemented with the provided data")
    }
    opt <- out[index,]
    
    do.call(c, list(as.list(out), opt=list(opt))) 
}

#====================================================================
#====================================================================
summary.SFI <- function(object,...)
{
    args0 <- list(...)
    if(!inherits(object, "SFI")) stop("The provided object is not from the class 'SFI'")
    
    tst <- object$tst
    y <- object$y
    if(object$method == "GBLUP"){
      df <- NA; lambda <- NA
    }else{
      df <- apply(object$df,2,mean)
      lambda <- apply(object$lambda,2,mean)
    }
    yHat <- as.matrix(fitted.SFI(object))
    
    correlation <- suppressWarnings(drop(stats::cor(y[tst],yHat,use="pairwise.complete.obs")))
    MSE <- suppressWarnings(apply((y[tst]-yHat)^2,2,sum,na.rm=TRUE)/length(tst))
    accuracy <- correlation/sqrt(object$h2)
    
    out <- data.frame(correlation=correlation, accuracy=correlation/sqrt(object$h2),
                MSE=MSE,df=df,lambda=lambda)

    # Detect maximum correlation
    index <- which.max(out$correlation)
    if(length(index)==0)
    {
      if(nrow(out)==1){
        index <- 1
      }else stop("'summary' method can not be implemented with the provided data")
    }
    opt <- out[index,]
  
    do.call(c, list(as.list(out), opt=list(opt))) 
}
