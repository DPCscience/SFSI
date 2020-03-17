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
saveBinary <- function(X,filename = paste0(tempdir(),"/file.bin"),
              row.names = TRUE, col.names = TRUE, size = 4, verbose = TRUE)
{
  if(!is.matrix(X)) stop("Object 'X' must be a matrix")
  if(storage.mode(X) != "double") storage.mode(X) <- "double"

  unlink(filename)
  if(row.names & !is.null(rownames(X))){
    rowNames <- rownames(X)
    if(any(nchar(rownames(X))>100)) stop("All rownames must be shorter than 100 characters long")
    nRowNames <- sum(nchar(rownames(X))) + nrow(X)
    if(nRowNames > .Machine$integer.max | nRowNames > 2^31 - 1) stop("Row names are too long")
  }else{
    rowNames <- ""
    nRowNames <- 0
  }
  if(col.names & !is.null(colnames(X))){
    colNames <- colnames(X)
    if(any(nchar(colnames(X))>100)) stop("All colnames must be shorter than 100 characters long")
    nColNames <- sum(nchar(colnames(X))) + ncol(X)
    if(nColNames > .Machine$integer.max | nColNames > 2^31 - 1) stop("Column names are too long")
  }else{
    colNames <- ""
    nColNames <- 0
  }
  out <- .Call('writeBinFile',filename,nrow(X),ncol(X),as.integer(size),X,
               as.integer(nRowNames),as.integer(nColNames),rowNames,colNames)
  if(verbose)
    cat("Saved file '",filename,"'\n nrows=",nrow(X)," ncols=",ncol(X)," size=",size," bytes\n")
}

#====================================================================
#====================================================================
readBinary <- function(filename = paste0(tempdir(),"/file.bin"),
                  indexRow = NULL, indexCol = NULL, verbose = TRUE)
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
  hasRowNames <- X[[4]]>0
  hasColNames <- X[[5]]>0
  rownamesX <- X[[6]]
  colnamesX <- X[[7]]
  if(nsetRow>0 & hasRowNames)  rownamesX <- rownamesX[indexRow]
  if(nsetCol>0 & hasColNames)  colnamesX <- colnamesX[indexCol]

  X <- X[[8]]
  if(hasRowNames) rownames(X) <- rownamesX
  if(hasColNames) colnames(X) <- colnamesX
  if(verbose)
    cat("Loaded file '",filename,"'\n nrows=",n," ncols=",p," size=",sizevar," bytes\n")

  return(X)
}

#====================================================================
# Derivative of the LogLikelihood in the ML
#====================================================================
dloglik <- function(lambda,n,Uty,UtX,d)
{
  # Trace of Hinv*G, Optimized
  dbar <- 1/(lambda*d+1)
  Tr_Hinv <- sum(dbar)
  Tr_Hinv_G <- (n-Tr_Hinv)/lambda

  qq1 <- t(Uty*dbar)%*%UtX
  qq2 <- solve(sweep(t(UtX),2L,dbar,FUN="*")%*%UtX)
  ytPy <- drop(sum(dbar*Uty^2)-qq1%*%qq2%*%t(qq1))

  # Check this, is suboptimal
  qq3 <- sum(Uty^2*dbar^2)
  qq4 <- drop(t(Uty*dbar^2)%*%UtX%*%qq2%*%t(qq1))
  qq5 <- drop(t(Uty*dbar)%*%UtX%*%qq2%*%t(UtX))*dbar
  qq6 <- sum(qq5^2)

  ytPPy <- drop(qq3-2*qq4+qq6)
  ytPGPy <- (ytPy-ytPPy)/lambda

  dL1 <- -0.5*Tr_Hinv_G + 0.5*n * ytPGPy/ytPy

  return(dL1)
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
# K=G; Z = indexK = subsetG = U = d = df = group.fill = group.shape = NULL; TST.col = "yellow"
# curve = FALSE; title = NULL; legend.pos = "topright"
# group = NULL; line.tick = 0.2; bg.col = "gray20"; line.col = "gray90"

plotNet <- function(fm, B, Z = NULL, K, indexK = NULL, subsetG = NULL,
           tst = NULL, U = NULL, d = NULL, group = NULL, group.fill = NULL,
           group.shape = NULL, group.size = NULL, df = NULL, title,
           curve = FALSE, TST.col = "yellow", bg.col = "gray20",
           line.col = "gray90", line.tick = 0.3, legend.pos="topright",
           set = c("Testing","Active","Non-active"))
{
  PC1 <- PC2 <- PC1_TRN <- PC1_TST <- PC2_TRN <- PC2_TST <- NULL
  legend.pos <- match.arg(legend.pos,
    choices=c("topright","bottomleft","bottomright","topleft","none"))

  if(class(fm) != "SFI") stop("Object 'fm' is not of the class 'SFI'")

  if(is.null(U) & is.null(d))
  {
    if(is.character(K)){
      K <- readBinary(K,indexRow=indexK,indexCol=indexK)
    }
    if(!is.double(K)) K <- apply(K,2,as.double)
    if(is.null(K))
      stop("A must be a positive semi definite matrix\n")
    if(!is.null(Z)) {
      if(!is.matrix(Z))  stop("Z must be a matrix\n")
      K <- Z %*% K %*% t(Z)
    }
    tmp <-  RSpectra::eigs_sym(K, 2)
    d <- tmp$values
    U <- tmp$vectors
    expvarPC <- 100*d/sum(diag(K))
  }else{
    if(is.null(U)){
      stop("You are providing the eigevalues, but not the eigenvectors")
    }else{
      if(is.null(d)){
        message("You are providing the eigenvectors, but not the eigenvalues\n",
                "No variance explained can be calculated")
        expvarPC <- NULL
      }else{
        if(nrow(U) == length(d)){
          expvarPC <- 100*(d^2)/sum(d^2)
        }else expvarPC <- NULL
      }
    }
  }
  tmp <- paste0(" (",sprintf('%.1f',expvarPC),"%)")
  if(length(tmp)<2) tmp <- NULL
  labelsPC <- paste0("PC",1:2,tmp[1:2])

  if(!is.null(tst)){
      if(any(!tst %in% fm$tst))
          stop("Any element in 'tst' vector is not contained in set 'fm$tst'")
  }else tst <- fm$tst

  justx <- ifelse(length(grep("left",legend.pos))>0,0,1)
  justy <- ifelse(length(grep("bottom",legend.pos))>0,0,1)
  legendPos <- "none"
  if(legend.pos != "none") legendPos <- c(abs(justx-0.01),abs(justy-0.01))

  theme0 <- theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.box.spacing = unit(0.4, "lines"),
    legend.background = element_rect(fill="gray95"),
    panel.background = element_rect(fill = bg.col),
    legend.justification = c(justx,justy),
    legend.position=legendPos,
    legend.key.height=unit(0.9,"line"),
    legend.key.width = unit(0.9, "lines")
  )

  if(missing(B)){
    if(is.null(df)) df <- summary(fm)$optCOR$df
    if(0 > df | df > range(fm$df)[2])
      stop("Parameter 'df' must be greater than zero and no greater than nTRN")
    B <- as.matrix(coef.SFI(fm,df=df))
  }else{
    stopifnot(is.matrix(B))
    df <- mean(apply(B,1,function(x)sum(abs(x)>0)))
  }

  flagGp <- !is.null(group)
  if(is.null(group)) group <- data.frame(group=rep(1,nrow(U)))
  gpName <- colnames(group)
  if(is.null(subsetG)) subsetG <- 1:nrow(U)

  if(!(class(set) == "character" & length(set) == 3))
   stop("Parameter 'set' must be a triplet of 'character' type")

  dat <- data.frame(id=1:nrow(U),set=set[3],group=group,U[,1:2])
  dat$set <- as.character(dat$set)

  # Testing and training (active) set
  dat$set[subsetG[tst]] <- set[1]
  index <- apply(B[fm$tst %in% tst,,drop=FALSE],2,function(x)any(abs(x)>0))
  dat$set[subsetG[fm$trn[index]]] <- set[2]
  dat$set[subsetG[fm$trn[!index]]] <- set[3]

  colnames(dat) <- c("id","set","group","PC1","PC2")

  dat$group <- factor(as.character(dat$group))
  dat$set <- factor(dat$set,levels=c(set))

  # Shape and color for the levels of group
  if(!flagGp) dat$group <- dat$set
  levelsGp <- levels(dat$group)

  if(is.null(group.fill)){
    if(flagGp){
      group.fill <- RColorBrewer::brewer.pal(8, "Set1")[1:7]
    }else group.fill <- c("#E69F00","#56B4E9","#999999")
    tmp <- ceiling(length(levelsGp)/length(group.fill))
    group.fill <- rep(group.fill,tmp)[1:length(levelsGp)]
  }
  if(is.null(group.shape)){
    if(flagGp){
      group.shape <- c(21,22,23,24,25)
    }else group.shape <- 21
    tmp <- ceiling(length(levelsGp)/length(group.shape))
    group.shape <- rep(group.shape,tmp)[1:length(levelsGp)]
  }
  if(is.null(group.size)){
    if(flagGp){
      group.size <- c(2.5,1,1)
    }else group.size <- c(2.5,1.5,1)
  }
  if(length(group.shape)!=length(levelsGp) | length(group.fill)!=length(levelsGp))
    stop("The number of elements in 'group.shape' and 'group.fill' must be of length ",length(levelsGp))

  if(missing(title)){
     title0 <- bquote(.(fm$name)*". Number of predictors="*.(round(df)))
     theme0 <- theme0 + theme(plot.title = element_text(hjust = 0.5))
  }else{
    title0 <- title
    if(is.null(title)){
      theme0 <- theme0 + theme(plot.title = element_blank())
    }else{
      theme0 <- theme0 + theme(plot.title = element_text(hjust = 0.5))
    }
  }

  indexTST <- which(dat$set == set[1])
  indexActive <- which(dat$set == set[2])
  indexNoActive <- which(!dat$set %in% set[1:2])

  tmp <- paste(paste(paste0("'",levelsGp,"'"),"=",group.shape),collapse=",")
  addShapeMan <- paste0("pt <- pt + scale_shape_manual(values=c(",tmp,"))")
  tmp <- paste(paste(paste0("'",levelsGp,"'"),"=",paste0("'",group.fill,"'")),collapse=",")
  addFillMan <- paste0("pt <- pt + scale_fill_manual(values=c(",tmp,"))")

  pt <- ggplot(dat,aes(x=PC1,y=PC2)) + labs(title=title0) +
         geom_point(data=dat[indexNoActive,],aes(shape=group,fill=group),
         color="black",size=group.size[3])

  for(i in 1:length(tst))
  {
    indexTRN <- which(abs(B[which(fm$tst == tst[i]),])>0)
    if(length(indexTRN)>0)
    {
      dat1 <- dat[subsetG[fm$trn],c("PC1","PC2")][indexTRN,]
      dat2 <- dat[subsetG[tst],c("PC1","PC2")][i,]
      colnames(dat1) <- paste0(colnames(dat1),"_TRN")
      colnames(dat2) <- paste0(colnames(dat2),"_TST")
      dat1 <- data.frame(dat2[rep(1,nrow(dat1)),],dat1)
      if(curve){
        pt <- pt + geom_curve(aes(x=PC1_TST,y=PC2_TST,xend=PC1_TRN,yend=PC2_TRN),
                      data=dat1,alpha=0.4,size=line.tick,color=line.col,curvature=0.4)
      }else{
        pt <- pt + geom_segment(aes(x=PC1_TST,y=PC2_TST,xend=PC1_TRN,yend=PC2_TRN),
                      data=dat1,alpha=0.4,size=line.tick,color=line.col)
      }
    }
  }

  pt <- pt  +
    geom_point(data=dat[dat$set==set[1],],aes(fill=group,shape=group),
               color=ifelse(flagGp,TST.col,"black"),size=group.size[1]) +
    geom_point(data=dat[dat$set==set[2],],aes(fill=group,shape=group),
               color="black",size=group.size[2]) +
    labs(shape=gpName,fill=gpName) + theme_bw() + theme0 +
    labs(x=labelsPC[1],y=labelsPC[2]) +
    guides(fill = guide_legend(override.aes = list(size = 2)))
  eval(parse(text=addShapeMan))
  eval(parse(text=addFillMan))

  if(!flagGp) pt <- pt + theme(legend.title = element_blank(),
                               legend.margin=margin(t=0,b=0.25,l=0.25,r=0.25,unit='line'))

  pt
}

#====================================================================
#====================================================================
plotPath <- function(fm, Z=NULL, K=NULL, indexK = NULL, tst=NULL, title=NULL, maxCor=0.85)
{
  k <- NULL
  flagKinship <- FALSE
  if(!class(fm) %in% c("SFI","SSI")) stop("Object 'fm' is not of the class 'SSI' or 'SFI'")

  theme0 <- theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.background = element_rect(fill="gray95"),
    legend.box.spacing = unit(0.4, "lines"),
    legend.key.height=unit(3,"line"),
    legend.key.width = unit(0.8, "lines")
  )

  if(class(fm) == "SFI")
  {
    if(!is.null(K)){
      flagKinship <- TRUE
      if(is.character(K)){
        K <- readBinary(K,indexRow=indexK,indexCol=indexK)
      }
      if(!is.double(K)) K <- apply(K,2,as.double)
      if(!is.null(Z)) K <- Z%*%K%*%t(Z)

      if(!is.matrix(K) | (length(K) != length(fm$y)^2))
        stop("Product Z %*% K %*% t(Z) must be a squared matrix with number of rows (and columns) equal to the number of elements in 'y'")
    }

    beta <- coef.SFI(fm)
    if(!is.null(tst)){
      if(any(!tst %in% fm$tst)) stop("Any element in 'tst' vector is not contained in set 'fm$tst'")
      indexTST <- which(fm$tst %in% tst)
    }else indexTST <- seq_along(fm$tst)
    beta <- beta[indexTST]
    lambda <- apply(fm$lambda,2,mean)
    df <- apply(fm$df,2,mean)

  }else{
    beta <- fm$beta
    lambda <- fm$lambda
    df <- fm$df
  }

  nDF <- length(df)
  if(min(lambda) < .Machine$double.eps*1000)  lambda[which.min(lambda)] <- min(lambda[lambda>0])/2
  if(nDF==1) stop("Coefficients path plot can not be generated for 'nLambda=1'")

  if(class(fm) == "SFI")
  {
    dat <- c()
    trim <- length(fm$trn)*length(beta) > 20000
    for(i in seq_along(beta))
    {
      b0 <- as.matrix(beta[[i]])
      if(trim){
        indexOK <- getIndexCorrelated(b0,maxCor)
      }else indexOK <- seq_along(fm$trn)

      tmp <- matrix(NA,nrow=1,ncol=length(indexOK))
      if(!is.null(K)) tmp <- K[fm$tst[indexTST[i]],fm$trn[indexOK],drop=FALSE]
      dimnames(tmp) <- list(fm$tst[indexTST[i]],fm$trn[indexOK])
      tmp <- reshape::melt(tmp)
      tmp <- tmp[rep(1:nrow(tmp),each=nDF),]

      df0 <- rep(df,length(indexOK))
      lambda0 <- rep(lambda,length(indexOK))
      b0 <- as.vector(b0[,indexOK])
      id <- factor(tmp$X1):factor(tmp$X2)
      dat <- rbind(dat,data.frame(df=df0,lambda=lambda0,beta=b0,k=tmp$value,id=id))
    }

  }else{
    id <- factor(rep(seq(ncol(beta)),each=nrow(beta)))
    dat <- data.frame(df=rep(df,ncol(beta)),lambda=rep(lambda,ncol(beta)),beta=as.vector(beta),id=id)
  }

  # Labels and breaks for the DF axis
  ax2 <- getSecondAxis(lambda,df)
  brks0 <- ax2$breaks
  labs0 <- ax2$labels

  title0 <- bquote("Coefficients path. "*.(fm$name))
  if(!is.null(title)) title0 <- title

  if(flagKinship)
  {
    pt <- ggplot(dat,aes(-log(lambda),beta,color=k,group=id)) +
      viridis::scale_color_viridis() + geom_line() + theme_bw() + theme0 +
      labs(title=title0,y=expression(beta),x=expression("-log("*lambda*")"))
  }else{
    pt <- ggplot(dat,aes(-log(lambda),beta,color=id,group=id))+
      geom_line() + theme_bw() + theme0 + theme(legend.position = "none") +
      labs(title=title0,y=expression(beta),x=expression("-log("*lambda*")"))
  }

  if(length(brks0)>3){
    pt <- pt + scale_x_continuous(sec.axis=sec_axis(~.+0,"Number of predictors",breaks=brks0,labels=labs0))
  }
  pt
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
  |  Sparse Family and Selection Index. Version 0.2.0 (Mar 16-2020)    |
  |====================================================================|
  ")
}
