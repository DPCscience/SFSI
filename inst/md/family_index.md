## Sparse Family Index (SFI)
In a Family Index the breeding value for each candidate of selection is estimated as a linear combination of the observed value of all the predictors (subjects in a training set). The contribution (regression coefficients) of all training subjects for each individual can be calculated simultaneously using the **BLUP (Best Linear Unbiased Predictor)** that relies in kinship relationship (either pedigree- or marker-based) between candidates of selection and training data. 

In contrast to the kinship-based BLUP, a Sparse Family Index (SFI) estimate the regression coefficients for each candidate with only a **subset** of the training subjects contributing to each individual's breeding value prediction. This predictor selection is achieved by imposing a penalization in the optimization function. The higher the value of the penalization parameter the smaller the number of predictors contributing to the prediction. The kinship-based BLUP appears as the un-penalized case of the SFI. 

Predictive ability of both kinship-based BLUP and SFI can be then compared using their **prediction accuracy** given by the correlation between observed and predicted values.

### Data
Data from CIMMYT’s Global Wheat Program. Lines were evaluated for grain yield (each entry corresponds to an average of two plot records) at four different environments; phenotypes (*wheat.Y* object) were centered and standardized to a unit variance within environment. Each of the lines were genotyped for 1279 diversity array technology (DArT) markers. At each marker two homozygous genotypes were possible and these were coded as 0/1. Marker genotypes are given in the object *wheat.X*. Finally a matrix *wheat.A* provides the pedigree relationships between lines computed from the pedigree records. Data is available for download in the R-package 'BGLR'.

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

**2. Comparing G-BLUP and non-sparse family index**

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
cbind(yHat_GBLUP,yHat_SFI)
```

Kinship-BLUP model can be fitted using function `GBLUP` from the 'SFSI' package
```r
fm <- GBLUP(G,y,h2)  
yHat_GBLUP2 <- predict(fm)$yHat  
head(cbind(yHat_GBLUP,yHat_SFI,yHat_GBLUP2))
```


**3. Comparing G-BLUP and non-sparse family index using cross-validation**

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

**4. Comparing G-BLUP and PFI for different values of the parameter lambda using cross-validation**

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

**5. Predicting values for a testing set using a training set**

```r
set.seed(123)
nTST <- 150               # Number of lines to predict
tst <- sample(1:n,nTST)   # Select lines to predict
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

Obtaining the value of `lambda` from a training data could be optimized by running repeatedly several cross-validations
by providing different values of the parameter `seed`

```r
# Repeated cross-validation in training data to get an optimal lambda
lambda <- c()
nRep <- 10   # Number of times to run the cross-validation
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

**6. Predicting values for a large testing set using parallel computing**

Analysis of a large number of individuals can be computational demanding. The options `nCores` and `subset` enable both parallel and distributed computing.
For parallel computing, option `mc.cores` allows to simultaneously run the program on several cores.
For distributed computing, `subset=c(j,nc)` divides the testing set into 'nc' chunks and run only the chunk 'j' separately. All the testing subsets can be separatelly run in a High Performance Computing (HPC) environment at different nodes. 

```r
tst <- sample(1:n,150)   # Select lines to predict
trn <- (1:n)[-tst]

nChunks <- 5    # Number of chunks in which the testing set will be divided into
j <- 1          # Subset to run at one node

# Run each of the subsets at different nodes and collect predicted values
for(j in 1:nChunks){
  fm <- SFI(G,y,h2,trn,tst,subset=c(j,nChunks),mc.cores=5)
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
plot(-log(lambda),correlation,type="l",lwd=2,col=2)
indexMax <- which.max(correlation)
abline(h=correlation[indexMax],lty=2,col=3)
abline(v=-log(lambda[indexMax]),lty=2,col=3)

# Remove files
unlink("out_subset_*.RData")
```

Results can be automatically saved in the 'working' directory at a provided path and prefix given in the `saveAt` parameter.

```r
prefix <- "testFolder/testSFI"      # Prefix (and path) that will be added to the output files name

# Run each of the subsets at different nodes
for(j in 1:nChunks){
  fm <- SFI(G,y,h2,trn,tst,subset=c(j,nChunks),saveAt=prefix,mc.cores=5)
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
yHat <- fitted(fm)
plot(fm)
plot(fm,G=G,PC=TRUE,df=10)     
```

The size and the number of output files might overflow disc memory, thus removing this files after use is advised
```r
unlink(paste0(prefix,"*.RData"))
unlink(paste0(prefix,"*.bin"))
```


