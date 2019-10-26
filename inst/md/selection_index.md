* [Back to main page](https://github.com/MarcooLopez/SFSI/blob/master/README.md)

## Sparse Selection Index (SSI)

### Data

Data will be simulated for *n* observations and *p* predictors. Both phenotypic values of the response variable *y* and predictors *X* are generated as the sum of some genotypic value plus some environmental deviation in such a way that there is a given correlation between the phenotypic and genotypic values (heritability). Also, some correlation exists between genotypic value of the response and that of all predictors (co-heritabilities), this value is equivalent to the squared root of the genetic correlation.

**1. Simulate data**

```r
set.seed(12345)
n <- 1000
p <- 1200

# Simulating response variable
h2y <- 0.3      
Uy = sqrt(h2y)*scale(rnorm(n))  # Genotypic value
Ey =  sqrt(1-h2y)*scale(rnorm(n))  # Environmental term
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
  Ux[,j] <- sqrt(h2x[j])*scale(a1 + a2)  # Genotypic value
  Ex[,j] <- sqrt(1-h2x[j])*scale(rnorm(n)) # Environmental term
  X[,j] <- scale(Ux[,j] + Ex[,j])
}
```

Genotypic covariances (between response and predictors) can be calculated from this simulated data by calculating the covariance between the simulated genotypic values; however in a real situation, genotypic values are not observed and covariances must be estimated from variance components using linear mixed models either using replicates or multivariate models considering kinship relationship among individuals
```r
# Genotypic covariance between response and predictors 
gencov <- drop(cov(Ux,Uy))
 
# Phenotypic covariance between response and predictors 
phencov <- drop(cov(X,y))

plot(phencov,gencov)

# Phenotypic covariance matrix among predictors
P <- var(X)
```

**2. Phenotypic vs Genotypic selection index**

```r
library(SFSI)

# Genotypic SS
fm1 <- SSI(P,gencov,method="CD")

# Phenotypic SS
fm2 <- SSI(P,phencov,method="CD")

# Regression coeficients for each value of lambda
B1 <- as.matrix(fm1$beta)
B2 <- as.matrix(fm2$beta)

# Fited values (selection indices)
yHat_Gen <- X %*% t(B1)
yHat_Phen <- X %*% t(B2)
```
The resulting indices are obtained such that the correlation between the index and the target is maximum (accuracy of selection). In this cases, the target of the phenotypic SS is the phenotype and for the genotypic SS is the genotype.

The accuracy of selection varies according to the penalization parameter lambda, thus an optimal value of lambda can be choosen such the accuracy of selection is maximum.

Again, the accuracy can be calculated using this simulated data but must be inferred from variance components in real data.

```r
corGen <- drop(cor(Uy,yHat_Gen))
corPhen <- drop(cor(Uy,yHat_Phen))

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
