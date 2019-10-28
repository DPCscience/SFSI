* [Back to main page](https://github.com/MarcooLopez/SFSI/blob/master/README.md)

## Sparse Selection Index (SSI)

### Data

Data will be simulated for *n* observations and *p* predictors. Both phenotypic values of the response variable *y* and predictors *X* are generated as the sum of some genotypic value plus some environmental deviation in such a way that there is a given correlation between the phenotypic and genotypic values (heritability). Also, some correlation exists between genotypic value of the response and that of all predictors (co-heritabilities), this value is equivalent to the squared root of the genetic correlation.

**1. Simulate data**

```r
set.seed(1234)
n <- 1500
p <- 1000

# Simulating response variable
h2y <- 0.3      
Uy = sqrt(h2y)*as.vector(scale(rnorm(n)))    # Genotypic value
Ey =  sqrt(1-h2y)*as.vector(scale(rnorm(n))) # Environmental term
y = scale(Uy + Ey)

# Co-heritabilities of the response with predictors
h2xy <- rbeta(p,3,8)

# Heritabilities of the predictors
h2x <- rbeta(p,8,8)

# Simulating predictor variables
Ux <- Ex <- x <- matrix(NA,ncol=p,nrow=n)

for(j in 1:p)
{
  a1 = sqrt(h2xy[j])*scale(Uy)
  a2 =  sqrt(1-h2xy[j])*scale(rnorm(n))
  Ux[,j] <- sqrt(h2x[j])*scale(a1 + a2)    # Genotypic value
  Ex[,j] <- sqrt(1-h2x[j])*scale(rnorm(n)) # Environmental term
  x[,j] <- scale(Ux[,j] + Ex[,j])
}
```

Genotypic covariances (between response and predictors) can be calculated from this simulated data by calculating the covariance between the simulated genotypic values; however in a real situation, genotypic values are not observed and covariances must be estimated from variance components using linear mixed models either using replicates or multivariate models considering kinship relationship among individuals
```r
# Genotypic covariance between response and predictors 
gencov <- as.vector(cov(Ux,Uy))
 
# Phenotypic covariance between response and predictors 
phencov <- as.vector(cov(x,y))

plot(phencov,gencov)

# Phenotypic covariance matrix among predictors
Px <- var(x)
```

**2. Phenotypic vs Genotypic selection index**

An index that uses phenotypic covariances will yield prediction that are the best in predicting phenotypic values; however this approach is not a good practice when the goal is selecting the best genotypes (judged by their genotypic values). In the later case, using genotypic covariances seems to be more appropiate.

Code below will calculate SS that use both phenotypic and genotypic covariances. 
```r
library(SFSI)

# Genotypic SS
fm1 <- SSI(Px,gencov,method="CD")

# Phenotypic SS
fm2 <- SSI(Px,phencov,method="CD")

# Regression coeficients for each value of lambda
B1 <- as.matrix(fm1$beta)
B2 <- as.matrix(fm2$beta)

# Fited values (selection indices)
GSI <- (x %*% t(B1))
PSI <- (x %*% t(B2))
```

The indices can be evaluated by their Mean Squared Error of prediction in genotypic value prediction
```r
# MSE of the indices
mseGSI <- apply((Uy-GSI)^2,2,sum)/n
msePSI <- apply((Uy-PSI)^2,2,sum)/n

dat <- rbind(
  data.frame(SI="GSI",MSE=mseGSI,df=fm1$df,lambda=fm1$lambda),
  data.frame(SI="PSI",MSE=msePSI,df=fm2$df,lambda=fm2$lambda)
)

library(ggplot2)
ggplot(dat[dat$df>1,],aes(-log(lambda),MSE,color=SI,group=SI)) + geom_line(size=0.8)
```
The resulting indices are obtained such that the correlation between the index and the target is maximum (accuracy of selection). In this cases, the target of the phenotypic SS is the phenotype and for the genotypic SS is the genotype.

**3. Phenotypic vs Genotypic selection index using cross-validation**

The accuracy of selection varies according to the penalization parameter lambda, thus an optimal value of lambda can be chosen such the accuracy of selection is maximum.

Again, the accuracy can be calculated using this simulated data but must be inferred from variance components in real data. The accuracy is equal to the product of the squared root of the heritability of the index times the genetic correlation between the index and the target 


```r
# Create folds to perform cross-validation
nFolds <- 5
seed <- 123
set.seed(seed)
folds <- rep(seq(1:nFolds), ceiling(n/nFolds))[1:n]
folds <- sample(folds)

# Objects to store results
accGSI <- accPSI <- c()
dfGSI <- dfPSI <- c()
lambdaGSI <- lambdaPSI <- c()

for(k in 1:nFolds)
{
  trn <- which(folds != k)
  tst <- which(folds == k)
  
  xTRN <- scale(x[trn,])
  yTRN <- scale(y[trn])
  
  UxTRN <- Ux[trn,]
  UyTRN <- Uy[trn]

  # Genotypic covariance between response and predictors 
  gencov <- drop(cov(UxTRN,UyTRN))

  # Phenotypic covariance between response and predictors 
  phencov <- drop(cov(xTRN,yTRN))

  # Phenotypic covariance matrix among predictors
  Px <- var(xTRN)

  fm1 <- SSI(Px,gencov,method="CD",tol=1E-4,maxIter=500)
  fm2 <- SSI(Px,phencov,method="CD",tol=1E-4,maxIter=500)
  
  # Retrieve data from 'df' and 'lambda'
  dfGSI <- cbind(dfGSI,fm1$df)
  dfPSI <- cbind(dfPSI,fm2$df)
  lambdaGSI <- cbind(lambdaGSI,fm1$lambda)
  lambdaPSI <- cbind(lambdaPSI,fm2$lambda)

  # Calculate the index (in testing set)
  xTST <- scale(x[tst,])
  GSI <- fitted(fm1,xTST)  
  PSI <- fitted(fm2,xTST)    
 
  # Accuracy of the indices (in testing set)
  UyTST <- Uy[tst]
  accGSI <- cbind(accGSI,drop(cor(UyTST,GSI)))
  accPSI <- cbind(accPSI,drop(cor(UyTST,PSI)))
  cat("----- Fold",k,". Done\n")
}

# Accuracy of the indices across folds
accGSIm <- apply(accGSI,1,mean)
accPSIm <- apply(accPSI,1,mean)

dfGSI <- apply(dfGSI,1,mean)
dfPSI <- apply(dfPSI,1,mean)
lambdaGSI <- apply(lambdaGSI,1,mean)
lambdaPSI <- apply(lambdaPSI,1,mean)

dat <- rbind(
  data.frame(SI="GSI",accuracy=accGSIm,df=dfGSI,lambda=lambdaGSI),
  data.frame(SI="PSI",accuracy=accPSIm,df=dfPSI,lambda=lambdaPSI)
)

# Plot the average accuracy (in testing set) across folds
ggplot(dat[dat$df>1,],aes(-log(lambda),accuracy,color=SI,group=SI)) + geom_line(size=0.8)

```

An optimal index can be obtained such as the accuracy is maximum
```r
dat <- rbind(
  data.frame(SI="GSI",accuracy=accGSI[which.max(accGSIm),]),
  data.frame(SI="PSI",accuracy=accPSI[which.max(accPSIm),])
)

ggplot(dat,aes(SI,accuracy)) + geom_bar(stat="identity",width=0.5)
```
