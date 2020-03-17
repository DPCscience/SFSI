* [Back to main page](https://github.com/MarcooLopez/SFSI/blob/master/README.md)

## Sparse Family Index (SFI)
In a Family Index the breeding value for each candidate of selection is estimated as a linear combination of the observed value of all the predictors (subjects in a training set). The contribution (regression coefficients) of all training subjects for each individual can be calculated simultaneously using the **BLUP (Best Linear Unbiased Predictor)** that relies in kinship relationship (either pedigree- or marker-based) between candidates of selection and training data.

In contrast to the kinship-based BLUP, a Sparse Family Index (SFI) estimate the regression coefficients for each candidate with only a **subset** of the training subjects contributing to each individual's breeding value prediction. This predictor selection is achieved by imposing a penalization in the optimization function. The higher the value of the penalization parameter the smaller the number of predictors contributing to the prediction. The kinship-based BLUP appears as the un-penalized case of the SFI.

Predictive ability of both kinship-based BLUP and SFI can be then compared using their **prediction accuracy** given by the correlation between observed and predicted values.

<div id="Outline" />

## Outline
  * [1. Data](#data)    
  * [2. Equivalence of G-BLUP and non-Sparse Family Index](#GBLUP&FI)
  * [3. Non-Sparse (G-BLUP) vs Sparse Family Index](#GBLUPvsSFI)
  * [4. Obtaining a value of lambda using cross-validation](#CV_SFI)
  * [5. Parallel computing for large datasets](#parallelizing)
  * [6. Working with binary files](#binaryFiles)

-------------------------------------------------------------------------------------------

<div id="data" />

### Data
Data from CIMMYT’s Global Wheat Program. Lines were evaluated for grain yield (each entry corresponds to an average of two plot records) at four different environments; phenotypes (*wheat.Y* object) were centered and standardized to a unit variance within environment. Each of the lines were genotyped for 1279 diversity array technology (DArT) markers. At each marker two homozygous genotypes were possible and these were coded as 0/1. Marker genotypes are given in the object *wheat.X*. Finally a matrix *wheat.A* provides the pedigree relationships between lines computed from the pedigree records. Data is available for download in the R-package 'BGLR' (Perez & de los Campos, 2014).

**1. Download and prepare data**

```r
library(SFSI)
data(wheat.E3)
data(wheat.G)

y <- tapply(X=Y$YLD,INDEX=as.character(Y$gid),FUN=mean)
y <- as.vector(y[rownames(G)])

# Calculate heritability
fm <- solveMixed(y,K=G)
varE <- fm$varE
varU <- fm$varU
h2 <- varU/(varU + varE)

```
[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="GBLUP&FI" />

**2. Equivalence of G-BLUP and non-Sparse Family Index**

A value of zero for the penalization parameter (`lambda=0`) yields an index whose regression coefficients are the same as those of the the genomic-BLUP model
```r
# G-BLUP  
G0 <- G
diag(G0) <- diag(G0) + (1-h2)/h2
B1 <- solve(G0)%*%G
uHat_GBLUP <- crossprod(B1,y-mean(y))        # Predicted values (in testing set)

# Non-sparse FI
fm <- SFI(y,K=G,h2=h2,lambda=0,mc.cores=4)  
B2 <- as.matrix(coef(fm))
uHat_SFI <- crossprod(B2,y-mean(y))          # Predicted values (in testing set)

# Compare regression coefficients
max(B1-B2)

# Compare results
cor(uHat_GBLUP,uHat_SFI)
plot(uHat_GBLUP,uHat_SFI)
head(out <- data.frame(uHat_GBLUP,uHat_SFI))

```
Slightly differences are due to the iterative feature of the SFI model. Adjusting parameters `tol` and `maxIter` can yield results
that perfectly match those of the G-BLUP model but default values provide sufficiently closed estimates.

***2.1. Fitting Kinship-based BLUP using 'solveMixed' function***

Predicted values can be directly retrieved from function `solveMixed` from the 'SFSI' package. The option `return.Hinv = TRUE` can be used to obtain the regression coefficients.

```r
fm <- solveMixed(y,K=G,return.Hinv=TRUE)  
B3 <- crossprod(G,fm$Hinv)
max(B1-B3)   # Regression coefficients comparison
uHat_GBLUP2 <- fm$u       # Predicted values (in testing set)
head(out <- data.frame(out,uHat_GBLUP2))
cor(out)

```

[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="GBLUPvsSFI" />

**3. Performance of G-BLUP (non-Sparse) and Sparse Family Index**

Predictive ability of the non-Sparse (or G-BLUP) will be compared with that of the Sparse Family Index using the correlation between
observed and predicted values in a cross-validation fashion for different values of the penalization parameter lambda.

Each value of `lambda` will yield a different family index.

```r
# Generate a grid of lambdas evenly spaced in logarithm scale starting from 1 to 0
nLambda <- 100     # Number of lambdas to generate
lambda <- exp(seq(log(1), log(.Machine$double.eps^0.5), length = nLambda))

# Create folds to perform cross-validation
nFolds <- 5
seed <- 123
set.seed(seed)
folds <- rep(seq(1:nFolds), ceiling(length(y)/nFolds))[seq_along(y)]
folds <- sample(folds)

corGBLUP <- rep(NA,ncol=nFolds)                # Object to store results of G-BLUP
corSFI <- matrix(NA,nrow=nLambda,ncol=nFolds)  # Object to store results of SFI

for(k in 1:nFolds)
{
  trn <- which(folds != k)
  tst <- which(folds == k)

  # G-BLUP
  yNA <- y
  yNA[tst] <- NA
  fm1 <- solveMixed(yNA,K=G)  
  uHat <- fm1$u[tst]           # Predicted values (in testing set)
  corGBLUP[k] <- cor(y[tst],uHat)

  # Penalized FI
  fm2 <- SFI(y,b=fm1$b,h2=fm1$h2,K=G,trn=trn,tst=tst,lambda=lambda,verbose=FALSE)  
  uHat <- fitted(fm2)           # Predicted values (in testing set)
  corSFI[,k] <- drop(cor(y[tst],uHat))
  cat("  -- Done fold=",k,"\n")
}

# Average across folds
corGBLUP <- mean(corGBLUP)
corSFI <- rowMeans(corSFI,na.rm=TRUE)   

# Plot results
plot(-log(lambda),corSFI,type="l",lwd=2,col=2,ylab="correlation")
abline(h=corGBLUP,lwd=2,col=3)
legend("bottomright",legend=c("SFI","GBLUP"),col=c(2,3),pch=20)

# Gain (in %) in accuracy of the maximum SFI
index <- which.max(corSFI)
100*(corSFI[index] - corGBLUP)/corGBLUP
```

***3.1 Comparison of performance using 'SFI_CV' function***

The above cross-validation can be done using the 'SFI_CV' function which divide data into a number of folds provided in `nFolds` parameter using the specified `seed` parameter.

```r
# Run the SFI with the generated grid of lambdas
fm1 <- SFI_CV(y,K=G,lambda=lambda,nFolds=nFolds,seed=seed,name="SFI",mc.cores=4)

# Run the un-penalized FI (G-BLUP model)
fm2 <- SFI_CV(y,K=G,lambda=0,nFolds=nFolds,seed=seed,name="G-BLUP",mc.cores=4)

# Plot of the (average) correlation in testing set vs the penalization parameter lambda
plot(fm1,fm2,py='accuracy')

# Plot of the (average) MSE in testing set vs the penalization parameter lambda
plot(fm1,fm2,py='MSE')

# Maximum average correlation
avgCor <- colMeans(do.call(rbind,lapply(fm1,function(x)(x$accuracy))))
max(avgCor,na.rm=TRUE)

# Maximum average correlation obtained using 'summary' method
(out1 <- summary(fm1)$optCOR["mean",])

# Average correlation obtained with G-BLUP
(out2 <- summary(fm2)$optCOR["mean",])

# Relative gain over G-BLUP (percentage)
100*(out1[3]-out2[3])/out2[3]
```

A vector of lambdas is generated internally by the program when `lambda` is not provided. The number of lambdas is passed as `nLambda` parameter.

```r
# Run the SFI without passing lambda but how many
fm3 <- SFI_CV(y,K=G,nLambda=nLambda,seed=seed,name="SFI",mc.cores=4)

# Plot of the (average) correlation in testing set vs the penalization parameter lambda
plot(fm3,fm2)
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/Figure_acc_vs_lambda_SFI.png" width="390">
</p>

[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="CV_SFI" />

**4. Predicting values for a testing set using a training set**

Breeding values can be predicted for individuals in a testing set using information from a training data.
An optimal value of the parameter `lambda` can be obtained by cross-validating in training data and then use it to
predict testing individuals.
```r
set.seed(1234)
pTST <- 0.3               # Proportion of lines to predict
tst <- sample(seq_along(y),ceiling(length(y)*pTST))   # Select lines to predict
trn <- (seq_along(y))[-tst]

# Cross-validation in training data to get a value of lambda
fm1 <- SFI_CV(y,K=G,trn.CV=trn,nFolds=5,mc.cores=4,seed=123,name="CV5")
lambda0 <- summary(fm1)$optCOR['mean','lambda']

# Predict testing data using lambda obtained from cross-validation
yNA <- y
yNA[tst] <- NA
fm <- SFI(yNA,K=G,trn=trn,tst=tst,lambda=lambda0,mc.cores=4)

# Correlation between predicted and observed values (in testing set)
plot(y[tst],fitted(fm))
cor(y[tst],fitted(fm))

# Correlation (in testing set) obtained with un-penalized FI (G-BLUP)
fm2 <- solveMixed(yNA,K=G)
cor(y[tst],fm2$u[tst])
```

An estimate of the parameter `lambda` could be also obtained by running repeatedly several cross-validations
and then averaging across repetitions. This task can be carried out either by setting `nCV` parameter
to a desired one or by providing a vector of different values through the parameter `seed`, in the latter case the number of
cross-validations to run are the length of the vector.

```r
fm2 <- SFI_CV(y,K=G,trn.CV=trn,nFolds=5,mc.cores=4,nCV=5,name="CV5x5")

# Lambda obtained across all fold/partitions
lambda0 <- summary(fm2)$optCOR['mean','lambda']             

fm <- SFI(yNA,K=G,trn=trn,tst=tst,lambda=lambda0)
cor(y[tst],fitted(fm))
```

Another cross-validation is the leave-one-out CV in which each individual is predicted using the remaining `n-1` individuals.
This CV is performed when setting `nFolds='n'`. This cross-validation might take long since it involves fitting `n-1` models

```r
fm3 <- SFI_CV(y,K=G,trn.CV=trn,nFolds='n',mc.cores=4,name="LOO")
lambda0 <- summary(fm3)$optCOR['mean','lambda']

fm <- SFI(yNA,K=G,trn=trn,tst=tst,lambda=lambda0)
cor(y[tst],fitted(fm))

# Comparison of the profile of each CV
plot(fm1,fm2,fm3,py="MSE")
plot(fm1,fm2,fm3)
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/CV_comparison_MSE_SFI.png" width="390">
</p>


***4.1 Network plot of individuals in testing and training***

The SFI gives for each individuals being predicted (testing set), an individualized set consisting of a reduced number of individuals (from training set) that contributes to its breeding value prediction. Function `plotNet` gives the representation of the conections (using lines) between testing and individuals in training that resulted with non-zero regression coefficient in the optimal SFI.

```r
plotNet(fm,K=G,bg.col="white",line.col="gray25")
       
# Passing a matrix of coefficients
B <- as.matrix(coef(fm))
tst0 <- fm$tst[1:15]   # A subset of the testing set
plotNet(fm,B=B,K=G,tst=tst0,curve=TRUE,group.size=c(3.5,1.5,1))
       
# Using Spectral Value Decomposition and grouping
EVD <- eigen(G)
gp <- data.frame(group=kmeans(EVD$vectors[,1:2],centers=5)$cluster)
plotNet(fm,tst=tst0,curve=TRUE,group=gp,U=EVD$vectors,d=EVD$values)
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/Network_plot_SFI.png" width="390">
</p>

[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="parallelizing" />

**5. Parallelizing computation of large number of testing individuals**

Analysis of a large number of individuals can be computationally demanding. The options `mc.cores` and `subset` enable both parallel and distributed computing.
For parallel computing, option `mc.cores` allows to simultaneously run the program on several cores.
For distributed computing, `subset=c(j,nc)` divides the testing set into 'nc' chunks and run only the chunk 'j' separately. All the testing subsets can be separatelly run in a High Performance Computing (HPC) environment at different nodes.

```r
set.seed(123)
pTST <- 0.2               # Proportion of lines to predict
tst <- sample(seq_along(y),ceiling(length(y)*pTST))   # Select lines to predict
trn <- (seq_along(y))[-tst]

nCores <- 4     # Number of cores in which the analysis will be run into
nChunks <- 5    # Number of chunks in which the testing set will be divided into
j <- 1          # Subset to run at one node

# Run each of the subsets at different nodes and collect predicted values
# for demonstration purposes all subsets will be run in a 'for' loop
for(j in 1:nChunks){
  fm <- SFI(y,K=G,trn=trn,tst=tst,subset=c(j,nChunks),mc.cores=nCores)
  uHat <- fitted(fm)
  df <- fm$df
  lambda <- fm$lambda
  yTST <- fm$y[fm$tst]
  save(uHat,yTST,df,lambda,file=paste0("out_subset_",j,".RData"))
}
```

Results can be collected after completion of all subsets
```r
nChunks <- 5
out <- vector('list',4)
for(j in 1:nChunks){
  load(paste0("out_subset_",j,".RData"))
  out$uHat <- rbind(out$uHat,uHat)
  out$yTST <- c(out$yTST,yTST)
  out$df <- rbind(out$df,df)
  out$lambda <- rbind(out$lambda,lambda)
}

# Accuracy in testing set
correlation <- drop(cor(out$yTST,out$uHat))
lambda <- apply(out$lambda,2,mean)
df <- apply(out$df,2,mean)
dat = data.frame(df=df,negLogLambda=-log(lambda),correlation)
plot(correlation~negLogLambda,data=dat[dat$df>1,],type="l",lwd=2,col=2)
indexMax <- which.max(correlation)
abline(h=correlation[indexMax],lty=2,col=3)
abline(v=-log(lambda[indexMax]),lty=2,col=3)

# Remove files
unlink("out_subset_*.RData")
```

Results can be automatically saved in the 'working' directory at a provided path and prefix given in the `saveAt` parameter.

```r
prefix <- "testSFI"      # Prefix (and path) that will be added to the output files name

# Run each of the subsets at different nodes
for(j in 1:nChunks){
  fm <- SFI(y,K=G,trn=trn,tst=tst,subset=c(j,nChunks),saveAt=prefix,mc.cores=nCores)
}
```

Providing `saveAt` parameter will generate `.RData` files where outputs are saved. Regression coefficients
are separatelly saved as binary (`*.bin`) files. Results of all chunks can be gathered after completion using function `collect`.

```r
fm <- collect(prefix)
```

Object `fm` does not contain the regression coefficients in memory but a path where they are storaged in disc. Methods `coef`, `summary`, `fitted`, and `plot` will read these files every time they are called

```r
df0 <- summary(fm)$optCOR$df  # Optimal value of degrees of freedom
uHat <- fitted(fm)
beta <- coef(fm)           # Coefficients for all values of lambda
beta <- coef(fm,df=df0)    # Coefficients for the optimum lambda
plot(fm)
plotNet(fm,K=G)            # Network plot
plotNet(fm,K=G,df=10)      # Network plot for a given 'df'

gp <- kmeans(eigen(G)$vectors[,1:2],centers=4)$cluster
gp <- data.frame(group = gp)
plotNet(fm,K=G,df=10,group=gp,TST.col="yellow")
plotPath(fm)               # Path plot   
plotPath(fm,K=G)           # Path plot using kinship  
```

The size and the number of output files might overflow disc memory, thus removing this files after use is advised
```r
unlink(paste0(prefix,"*.RData"))
unlink(paste0(prefix,"*.bin"))
```
[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="binaryFiles" />

**6. Working with binary files**

Large matrices occupy a big space on disk and hence take long time to be loaded into memory in R. Genomic matrix can be saved as binary
file as single-precision *float* type (4 bytes, 32-bits) of space and thus reducing to half the space usage relative to the double-precision *double* type (8 bytes, 64-bits). This can be done using function `saveBinary` specifying either 4 or 8 bytes in parameter `size`.
```r
# Save G matrix as binary file
saveBinary(G,"G_matrix_32bits.bin",size=4)   # as single-precision
saveBinary(G,"G_matrix_64bits.bin",size=8)   # as double-precision

# Size of files (in Megabytes)
file.size("G_matrix_32bits.bin")/(1024^2)
file.size("G_matrix_64bits.bin")/(1024^2)
```

Saving the matrix at 4 bytes occupies less space at the cost of reducing the precision relative to 8 bytes that *numeric* variables occupy in *R*.
```r
# Read G matrix from a previously saved binary file
G2 <- readBinary("G_matrix_32bits.bin")   
G3 <- readBinary("G_matrix_64bits.bin")

# Numeric precision
sum(abs(G-G2))  # Loss of precision relative to the original matrix
sum(abs(G-G3))  # No loss of precision
```

***6.1. Passing a binary file name instead of a matrix***

Parameter `G` can be the name of a binary file containing a genomic matrix which will be read internally. In addition, specific rows
and columns can be read using parameter `indexG` which should match to the individuals whose scores are passed in `y`.
```r
fm <- SFI(y,K="G_matrix_32bits.bin",trn=trn,tst=tst)
fm <- SFI_CV(y,K="G_matrix_32bits.bin")

# Selecting specific individuals to work with
y2 <- y[1:500]             # Work only with first 500 individuals
tst <- sample(1:length(y2),floor(length(y2)*0.2))   # Select 20% of lines to predict
trn <- (1:length(y2))[-tst]
fm <- SFI(y2,K="G_matrix_32bits.bin",trn=trn,tst=tst,indexK=1:500)
fm <- SFI_CV(y2,K="G_matrix_32bits.bin",indexK=1:500)
```

[Back to Outline](#Outline)
