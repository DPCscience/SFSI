## Penalized Family Index
In a Family Index the breeding value for each candidate of selection is estimated as a linear combination of the observed value of all the predictors (subjects in a training set). The contribution (regression coefficients) of all training subjects for each individual can be calculated simultaneously using the **BLUP (Best Linear Unbiased Predictor)** that relies in kinship relationship (either pedigree- or marker-based) between candidates of selection and training data. 

In contrast to the kinship-based BLUP, the penalized Family Index (PFI) estimate the regression coefficients for each candidate with only a **subset** of the training subjects contributing to each individual's breeding value prediction. The higher the value of the penalization parameter the smaller the number of predictors contributing to the prediction. The kinship-based BLUP appears as the un-penalized case of the PFI. 

Predictive ability of both kinship-based BLUP and PFI can be then compared using their prediction accuracy given by the correlation between observed and predicted values.

### 1. Data
Data from CIMMYT’s Global Wheat Program. Lines were evaluated for grain yield (each entry corresponds to an average of two plot records) at four different environments; phenotypes (*wheat.Y* object) were centered and standardized to a unit variance within environment. Each of the lines were genotyped for 1279 diversity array technology (DArT) markers. At each marker two homozygous genotypes were possible and these were coded as 0/1. Marker genotypes are given in the object *wheat.X*. Finally a matrix *wheat.A* provides the pedigree relationships between lines computed from the pedigree records. Data is available for download in the R-package 'BGLR'.

**Download data**
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
```

### 2. Fitting G-BLUP and penalized family index
```r
# Calculating heritability using rrBLUP package
# install.packages("rrBLUP")  # If not installed
library(rrBLUP)
In <- diag(n)
fm <- mixed.solve(y=y,Z=In,K=G)
varE <- fm$Ve
varU <- fm$Vu
h2 <- varU/(varU + varE)

# Creating folds to perform cross-validation
nFolds=3
set.seed(123)
folds <- rep(seq(1:nFolds), ceiling(n/nFolds))[1:n]
folds <- sample(folds)

# Calculating G-BLUP
for(k in 1:nFolds)
{
  trn <- which(folds != k)
  tst <- which(folds == k)
  fm <- PFI(G,y,h2,trn,tst,verbose=TRUE)
  fm <- mixed.solve(y=y,Z=In,K=G)
}


for(k in 1:nFolds)
{
  trn <- which(folds != k)
  tst <- which(folds == k)
  fm <- PFI(G,y,h2,trn,tst,verbose=TRUE)
}


```
