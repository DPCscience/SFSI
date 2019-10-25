* [Back to main page](https://github.com/MarcooLopez/SFSI/blob/master/README.md)

## Sparse Selection Index (SSI)

### Data

**1. Simulate data**

```r
set.seed(12345)
n <- 1000
p <- 750

# Simulating response variable
h2y <- 0.4
Uy = sqrt(h2y)*scale(rnorm(n))
Ey =  sqrt(1-h2y)*scale(rnorm(n))
y = scale(Uy + Ey)

# Co-heritabilities of the response with predictors
h2xy <- rbeta(p,3,8)

# Heritabilities of the predictors
h2x <- rbeta(p,8,3)

# Simulating predictor variables
Ux <- Ex <- X <- matrix(NA,ncol=p,nrow=n)

for(j in 1:p)
{
  a1 = sqrt(h2xy[j])*scale(Uy)
  a2 =  sqrt(1-h2xy[j])*scale(rnorm(n))
  Ux[,j] <- sqrt(h2x[j])*scale(a1 + a2)
  Ex[,j] <- sqrt(1-h2x[j])*scale(rnorm(n))
  
  X[,j] <- scale(Ux[,j] + Ex[,j])
}
```

Phenotypic and genetic correlations can be calculated from this simulated data; however for real data, genotypic covariance
must be estimated thorugh linear mixed models using replicates or using kinship relationship among individuals
```r
# Phenotypic covariance between response and predictors 
phencov <- drop(cov(X,y))

# Genotypic covariance between response and predictors 
gencov <- drop(cov(Ux,Uy))
  
plot(phencov,gencov)

# Phenotypic covariance matrix among predictors
P <- var(X)
```

**2. Phenotypic vs Genotypic selection index**

```r
library(SFSI)

fm1 <- SSI(P,gencov,method="CD")
fm2 <- SSI(P,phencov,method="CD")

yHatGen <- fitted(fm1,X)
yHatPhen <- fitted(fm2,X)

corPhen <- drop(cor(Uy,yHatPhen))
corGen <- drop(cor(Uy,yHatGen))

rg <- range(c(corPhen,corGen),na.rm=TRUE)

plot(-log(fm1$lambda),corGen,col=2,type="l",lwd=2,xlab=expression(-log(lambda)),ylab="accuracy",ylim=rg)
points(-log(fm2$lambda),corPhen,col=4,type="l",lwd=2)
legend("bottomleft",legend=c("GenSS","PhenSS"),col=c(2,4),pch=20)
```

**3. Phenotypic vs Genotypic selection index using cross-validation**

```r

# Create folds to perform cross-validation
nFolds <- 10
seed <- 123
set.seed(seed)
folds <- rep(seq(1:nFolds), ceiling(n/nFolds))[1:n]
folds <- sample(folds)

for(k in 1:nFolds)
{
  trn <- which(folds != k)
  tst <- which(folds == k)

  # Phenotypic covariance between response and predictors 
  phencov <- drop(cov(X[trn,],y[trn]))

  # Genotypic covariance between response and predictors 
  gencov <- drop(cov(Ux[trn,],Uy[trn]))
  
  # Phenotypic covariance matrix among predictors
  P <- var(X[trn,])

  fm1 <- SSI(P,gencov,method="CD")
  fm2 <- SSI(P,phencov,method="CD")

  yHatGen <- fitted(fm1,X[trn,])
  yHatPhen <- fitted(fm2,X[trn,])

  corPhen <- drop(cor(Uy,yHatPhen))
  corGen <- drop(cor(Uy,yHatGen))
}
rg <- range(c(corPhen,corGen),na.rm=TRUE)

plot(fm1$df,corGen,col=2,type="l",lwd=2,xlab="number of predictors",ylab="accuracy",ylim=rg)
points(fm2$df,corPhen,col=4,type="l",lwd=2)
legend("bottomleft",legend=c("GenSS","PhenSS"),col=c(2,4),pch=20)

```
