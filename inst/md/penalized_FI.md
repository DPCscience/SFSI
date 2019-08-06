## Penalized Family Index (PFI)
In a Family Index the breeding value for each candidate of selection is estimated as a linear combination of the observed value of all the predictors (subjects in a training set). The contribution (regression coefficients) of all training subjects for each individual can be calculated simultaneously using the **BLUP (Best Linear Unbiased Predictor)** that relies in kinship relationship (either pedigree- or marker-based) between candidates of selection and training data. 

In contrast to the kinship-based BLUP, the penalized Family Index (PFI) estimate the regression coefficients for each candidate with only a **subset** of the training subjects contributing to each individual's breeding value prediction. The higher the value of the penalization parameter the smaller the number of predictors contributing to the prediction. The kinship-based BLUP appears as the un-penalized case of the PFI. 

Predictive ability of both kinship-based BLUP and PFI can be then compared using their prediction accuracy given by the correlation between observed and predicted values.

### Data
Data from CIMMYT’s Global Wheat Program. Lines were evaluated for grain yield (each entry corresponds to an average of two plot records) at four different environments; phenotypes (*wheat.Y* object) were centered and standardized to a unit variance within environment. Each of the lines were genotyped for 1279 diversity array technology (DArT) markers. At each marker two homozygous genotypes were possible and these were coded as 0/1. Marker genotypes are given in the object *wheat.X*. Finally a matrix *wheat.A* provides the pedigree relationships between lines computed from the pedigree records. Data is available for download in the R-package 'BGLR'.

**1. Download and prepare data**
```r
# install.packages("BGLR")  # If not installed
library(BGLR)
data(wheat)
A <- wheat.A
X <- wheat.X
Y <- wheat.Y

# Select an environment
y <- as.vector(Y[,1])
n <- length(y)

# Calculating the genomic relationship matrix
G <- tcrossprod(scale(X))/ncol(X)

# Calculating heritability using 'rrBLUP' package
# install.packages("rrBLUP")  # If not installed
library(rrBLUP)
In <- diag(n)
fm <- mixed.solve(y=y,Z=In,K=G)
varE <- fm$Ve
varU <- fm$Vu
h2.0 <- varU/(varU + varE)

# Creating folds to perform cross-validation
nFolds=4
set.seed(123)
folds <- rep(seq(1:nFolds), ceiling(n/nFolds))[1:n]
folds <- sample(folds)
```

**2. Comparing G-BLUP and un-penalized family index using cross-validation**
```r
library(PFSI)
out <- matrix(NA,ncol=3,nrow=nFolds)    # Object to store results
colnames(out) <- c("rrBLUP","PFI1","PFI2")
h2 <- c()         # To save within-fold heritability
for(k in 1:nFolds)
{
  trn <- which(folds != k)
  tst <- which(folds == k)
  yNA <- y
  yNA[tst] <- NA
  
  # G-BLUP using 'rrBLUP' package
  fm <- mixed.solve(y=yNA,Z=In,K=G)
  out[k,'rrBLUP'] <- cor(fm$u[tst],y[tst])
  h2[k] <- fm$Vu/(fm$Vu + fm$Ve)
  
  # G-BLUP as un-penalized FI using within-fold heritability
  fm <- PFI(G,y,h2[k],trn,tst,lambda=0)
  out[k,'PFI1'] <- cor(predict(fm)$yHat,y[tst])
  
  # G-BLUP as un-penalized FI using using heritability calculated using complete data
  fm <- PFI(G,y,h2.0,trn,tst,lambda=0)
  out[k,'PFI2'] <- cor(predict(fm)$yHat,y[tst])
  cat("Done fold=",k,"\n")
}

# Comparing results
out
apply(out,2,mean)
```
  
The above cross-validation can be done using the 'PFI_CV' function with the same 'seed' and same number of folds (`seed=123` and `nFolds=4`)

```r
fm <- PFI_CV(G,y,h2.0,lambda=0,nFolds=4,seed=123)

# Comparing with previous results (that used heritability calculated using complete data)
cbind(fm$correlation,out[,'PFI2'])
```

**3. Comparing G-BLUP and PFI for different values of the parameter lambda**
```r

# Generating a grid of 100 lambdas evenly spaced in logarithm scale starting from 1 to 0
nLambda <- 100
lambda <- exp(seq(log(1), log(1e-05), length = nLambda))
lambda[nLambda] <- 0

# Running the PFI with the generated lambdas
fm <- PFI_CV(G,y,h2.0,lambda=lambda,nFolds=4,seed=123,name="PFI")

# Running the un-penalized FI
fm2 <- PFI_CV(G,y,h2.0,lambda=0,nFolds=4,seed=123,name="G-BLUP")

# Plot of the (average) correlation in testing set vs the penalization parameter lambda
plot(fm1,fm2,px='lambda')

# Plot of the (average) correlation in testing set vs the average number of predictors (in training set)
plot(fm1,fm2,px='df')

# Getting the maximum correlation
avgCor <- apply(fm1$correlation,2,mean)
max(avgCor,na.rm=TRUE)

# Optimal penalization that gives the highest correlation
summary(fm1)[['PFI']][[1]][['peaks']]

# Relative gain over G-BLUP
summary(fm2)[['PFI']][[1]][['peaks']][1,'correlation']

In <- diag(n)
corGBLUP <- c()
h2 <- c()
for(k in 1:nFolds)
{
  yNA <- y
  tst <- which(folds == k)
  yNA[tst] <- NA
  fm <- mixed.solve(y=yNA,Z=In,K=G)
  corGBLUP[k] <- cor(fm$u[tst],y[tst])
  h2[k] <- fm$Vu/(fm$Vu + fm$Ve)
}

# Calculating G-BLUP as a un-penalized family index (lambda=0) using 'PFSI' package
library(PFSI)
corPFI <- c()
for(k in 1:nFolds)
{
  trn <- which(folds != k)
  tst <- which(folds == k) 
  fm <- PFI(G,y,h2[k],trn,tst,lambda=0)
  corPFI[k] <- cor(predict(fm)$yHat,y[tst])
}

# Comparing both results
cbind(corGBLUP,corPFI)
mean(corGBLUP);mean(corPFI)
```
