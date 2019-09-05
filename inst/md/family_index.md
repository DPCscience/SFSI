## Sparse Family Index (SFI)
In a Family Index the breeding value for each candidate of selection is estimated as a linear combination of the observed value of all the predictors (subjects in a training set). The contribution (regression coefficients) of all training subjects for each individual can be calculated simultaneously using the **BLUP (Best Linear Unbiased Predictor)** that relies in kinship relationship (either pedigree- or marker-based) between candidates of selection and training data. 

In contrast to the kinship-based BLUP, a Sparse Family Index (SFI) estimate the regression coefficients for each candidate with only a **subset** of the training subjects contributing to each individual's breeding value prediction. This predictor selection is achieved by imposing a penalization in the optimization function. The higher the value of the penalization parameter the smaller the number of predictors contributing to the prediction. The kinship-based BLUP appears as the un-penalized case of the SFI. 

Predictive ability of both kinship-based BLUP and SFI can be then compared using their **prediction accuracy** given by the correlation between observed and predicted values.

<div id="Outline" />

## Outline
  * [Data](#data)    
  * [Equivalence of G-BLUP and non-Sparse Family Index](#GBLUP&FI)
  * [non-Sparse (G-BLUP) vs Sparse Family Index](#GBLUPvsSFI)
  * [Predicting testing individuals using a trainig set](#predictSFI)
  * [Parallel computing for large datasets](#parallelizing)
  * [Working with binary files](#binaryFiles)
   
-------------------------------------------------------------------------------------------

<div id="data" />

### Data
Data from CIMMYT’s Global Wheat Program. Lines were evaluated for grain yield (each entry corresponds to an average of two plot records) at four different environments; phenotypes (*wheat.Y* object) were centered and standardized to a unit variance within environment. Each of the lines were genotyped for 1279 diversity array technology (DArT) markers. At each marker two homozygous genotypes were possible and these were coded as 0/1. Marker genotypes are given in the object *wheat.X*. Finally a matrix *wheat.A* provides the pedigree relationships between lines computed from the pedigree records. Data is available for download in the R-package 'BGLR' (Perez & de los Campos, 2014).

**1. Download and prepare data**

```r
data(wheat,package="BGLR")
A <- wheat.A
X <- wheat.X
Y <- wheat.Y

# Select an environment
y <- as.vector(Y[,1])
n <- length(y)

# Calculate the genomic relationship matrix
G <- tcrossprod(scale(X))/ncol(X)

# Calculate heritability using 'BGLR' package
fm <- BGLR::BGLR(y,ETA=list(list(K=G,model="RKHS")),nIter=12000,burnIn=5000)
varE <- fm$varE
varU <- fm$ETA[[1]]$varU
h2 <- varU/(varU + varE)

```
[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="GBLUP&FI" />

**2. Equivalence of G-BLUP and non-Sparse Family Index**

A value of zero for the penalization parameter (`lambda=0`) yields an index whose regression coefficients are the same as those of the the genomic-BLUP model
```r
library(SFSI)

# G-BLUP  
G0 <- G
diag(G0) <- diag(G0) + (1-h2)/h2
B <- solve(G0)%*%G
yHat_GBLUP <- crossprod(B,y-mean(y))         # Predicted values (in testing set)
  
# Non-sparse FI
fm <- SFI(G,y,h2,lambda=0,mc.cores=4,verbose=TRUE)  
yHat_SFI <- predict(fm)$yHat                 # Predicted values (in testing set)

# Compare regression coefficients
max(B-coef(fm,df=length(fm$training)))

# Compare results
cor(yHat_GBLUP,yHat_SFI)
plot(yHat_GBLUP,yHat_SFI)
head(cbind(yHat_GBLUP,yHat_SFI))
```
Slightly differences are due to the iterative feature of the SFI model. Adjusting parameters `tol` and `maxIter` can yield results
that perfectly match those of the G-BLUP model but default values provide sufficiently closed estimates. 

Kinship-BLUP model can be fitted using function `GBLUP` from the 'SFSI' package
```r
fm <- GBLUP(G,y,h2)  
yHat_GBLUP2 <- predict(fm)$yHat  
head(cbind(yHat_GBLUP,yHat_SFI,yHat_GBLUP2))
```

**2.1. Equivalence of G-BLUP and non-Sparse Family Index using cross-validation**

```r
# Create folds to perform cross-validation
nFolds <- 5
seed <- 123
set.seed(seed)
folds <- rep(seq(1:nFolds), ceiling(n/nFolds))[1:n]
folds <- sample(folds)

out <- matrix(NA,ncol=2,nrow=nFolds)    # Object to store results
colnames(out) <- c("GBLUP","SFI")
for(k in 1:nFolds)
{
  trn <- which(folds != k)
  tst <- which(folds == k)
  
  # G-BLUP 
  fm <- GBLUP(G,y,h2,trn,tst)  
  yHat <- predict(fm)$yHat  
  out[k,'GBLUP'] <- cor(yHat,y[tst])
  
  # G-BLUP as un-penalized FI
  fm <- SFI(G,y,h2,trn,tst,lambda=0)  
  yHat <- predict(fm)$yHat               # Predicted values (in testing set)
  out[k,'SFI'] <- cor(yHat,y[tst])
  cat("  -- Done fold=",k,"\n")
}

# Compare results
out
colMeans(out)   # Average across folds
```
  
The above cross-validation can be done using the 'PFI_CV' function with the same `seed` and `nFolds` parameters

```r
fm <- SFI_CV(G,y,h2,lambda=0,nFolds=nFolds,seed=seed,mc.cores=4)

# Compare with previous results (that used heritability calculated using complete data)
cbind(out,fm$correlation)
```
[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="GBLUPvsSFI" />

**3. Performance of G-BLUP (non-Sparse) and Sparse Family Index**

Predictive ability of the non-Sparse (or G-BLUP) will be compared with that of the Sparse Family Index using the correlation between
observed and predicted values in a cross-validation fashion for different values of the penalization parameter lambda.
```r
# Generate a grid of lambdas evenly spaced in logarithm scale starting from 1 to 0
nLambda <- 100     # Number of lambdas to generate
lambda <- exp(seq(log(1), log(1e-05), length = nLambda))
lambda[nLambda] <- 0

# Run the SFI with the generated grid of lambdas
fm1 <- SFI_CV(G,y,h2,lambda=lambda,nFolds=nFolds,seed=seed,name="SFI",mc.cores=4)

# Run the un-penalized FI (G-BLUP model)
fm2 <- SFI_CV(G,y,h2,lambda=0,nFolds=nFolds,seed=seed,name="G-BLUP",mc.cores=4)

# Plot of the (average) correlation in testing set vs the penalization parameter lambda
plot(fm1,fm2,py='correlation')

# Plot of the (average) MSE in testing set vs the penalization parameter lambda
plot(fm1,fm2,py='MSE')

# Maximum average correlation
avgCor <- colMeans(fm1$correlation)
max(avgCor,na.rm=TRUE)

# Maximum average correlation obtained using 'summary' method 
out1 <- summary(fm1)[['SFI']][[1]][['max']]
out1

# Maximum average correlation for G-BLUP
out2 <- summary(fm2)[['SFI']][[1]][['max']]
out2

# Relative gain over G-BLUP (percentage)
100*(out1[1]-out2[1])/out1[1]
```

The same comparison between G-BLUP and PFI can be done without passing a vector of lambdas since they are generated internally
by the program when `lambda` is not provided

```r
# Run the PFI. Lambdas will be generated automatically
fm1 <- SFI_CV(G,y,h2,nFolds=nFolds,seed=seed,mc.cores=4)

# Run the un-penalized FI. Using option method='G-BLUP'
fm2 <- SFI_CV(G,y,h2,method="GBLUP",nFolds=nFolds,seed=seed)

# Plot of the (average) correlation in testing set vs the average number of predictors (in training set)
plot(fm1,fm2)

# Relative gain over G-BLUP (percentage)
summary(fm1,fm2)[['SFI']][[1]][['gain']]
```
[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="predictSFI" />

**4. Predicting values for a testing set using a training set**

Breeding values can be predicted for individuals in a testing set using information from a training data.
An optimal value of the parameter `lambda` can be obtained by cross-validating in training dat and then use it to
predict testing individuals.
```r
set.seed(123)
pTST <- 0.2               # Proportion of lines to predict
tst <- sample(1:n,floor(n*pTST))   # Select lines to predict
trn <- (1:n)[-tst]

# Cross-validation in training data to get a value of lambda
fm <- SFI_CV(G,y,h2,training=trn,nFolds=3,mc.cores=4,seed=123)
lambda <- summary(fm)[['SFI']][[1]][['max']][1,'lambda']

# Predict testing data using lambda obtained from cross-validation
yNA <- y
yNA[tst] <- NA
fm <- SFI(G,yNA,h2,trn,tst,lambda=lambda,mc.cores=4)

# Correlation between predicted and observed values (in testing set)
plot(predict(fm)$yHat,y[tst])
cor(predict(fm)$yHat,y[tst])

# Correlation (in testing set) obtained with un-penalized FI (G-BLUP)
GBLUP(G,y,h2,trn,tst)$correlation
```

A better estimate of the parameter `lambda` could be obtained by running repeatedly several cross-validations
by providing different values of the parameter `seed` and then averaging across repetitions.

```r
# Repeated cross-validation in training data to get an optimal lambda
nRep <- 10      # Number of times to run the cross-validation
lambda <- c()   # Vector to store the values of lambda
for(j in 1:nRep)
{
   fm <- SFI_CV(G,y,h2,training=trn,nFolds=3,mc.cores=4,seed=j*500)
   lambda[j] <- summary(fm)[[1]][[1]][['max']][1,'lambda']
   cat("  -- Done repetition=",j,"\n")
}

# Obtain an optimal lambda by averaging the ones obtained by cross-validation
lambda0 <- mean(lambda)                     # Aritmethic mean or
lambda0 <- prod(lambda)^(1/length(lambda))  # Geometric mean

fm <- SFI(G,yNA,h2,trn,tst,lambda=lambda0)
cor(predict(fm)$yHat,y[tst])
```
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
tst <- sample(1:n,floor(n*pTST))   # Select lines to predict
trn <- (1:n)[-tst]

nCores <- 4     # Number of cores in which the analysis will be run into
nChunks <- 5    # Number of chunks in which the testing set will be divided into
j <- 1          # Subset to run at one node

# Run each of the subsets at different nodes and collect predicted values
# for demonstration purposes all subsets will be run in a 'for' loop
for(j in 1:nChunks){
  fm <- SFI(G,y,h2,trn,tst,subset=c(j,nChunks),mc.cores=nCores)
  yHat <- fitted(fm)
  df <- fm$df
  lambda <- fm$lambda
  yTST <- fm$y[fm$testing] 
  save(yHat,yTST,df,lambda,file=paste0("out_subset_",j,".RData"))
}
```

Results can be collected after completion of all subsets
```r
nChunks <- 5 
out <- vector('list',4)
for(j in 1:nChunks){
  load(paste0("out_subset_",j,".RData"))
  out$yHat <- rbind(out$yHat,yHat)
  out$yTST <- c(out$yTST,yTST)
  out$df <- rbind(out$df,df)
  out$lambda <- rbind(out$lambda,lambda)
}

# Accuracy in testing set
correlation <- drop(cor(out$yTST,out$yHat))
lambda <- apply(out$lambda,2,mean)
dat = data.frame(df=apply(out$df,2,mean),logLambda=-log(lambda),correlation)
plot(correlation~logLambda,data=dat[dat$df>1,],type="l",lwd=2,col=2)
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
  fm <- SFI(G,y,h2,trn,tst,subset=c(j,nChunks),saveAt=prefix,mc.cores=nCores)
}
```

Providing `saveAt` parameter will generate `.RData` files where outputs are saved. Regression coefficients
are separatelly saved as binary (`*.bin`) files. Results of all chunks can be gathered after completion using function `collect`. 

```r
fm <- collect(prefix)
```

Object `fm` does not contain the regression coefficients in memory but a path where they are storaged in disc. Methods `coef`, `summary`, `predict`, `fitted`, and `plot` will read these files every time they are called

```r
summary(fm)
df <- summary(fm)[[1]][['max']][1,'df']  # Optimal value of degrees of freedom
yHat <- fitted(fm)
beta <- coef(fm)       # Coefficients for all values of lambda
beta <- coef(fm,df=df) # Coefficients for the optimum lambda
plot(fm)
plot(fm,G=G)     # Coefficients path   
plot(fm,G=G,PC=TRUE,df=10)     
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

**6.1. Passing a binary file name instead of a matrix**

Parameter `G` can be the name of a binary file containing a genomic matrix which will be read internally. In addition, specific rows
and columns can be read using parameter `indexG` which should match to the individuals whose scores are passed in `y`.
```r
fm <- SFI("G_matrix_32bits.bin",y,h2,trn,tst)
summary(fm)

# Selecting specific individuals to work with
y2 <- y[1:500]             # Work only with first 500 individuals
tst <- sample(1:length(y2),floor(length(y2)*0.2))   # Select 20% of lines to predict
trn <- (1:length(y2))[-tst]
fm <- SFI("G_matrix_32bits.bin",y2,h2,trn,tst,indexG=1:500)
summary(fm)
```

[Back to Outline](#Outline)

