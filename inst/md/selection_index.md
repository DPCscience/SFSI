* [Back to main page](https://github.com/MarcooLopez/SFSI/blob/master/README.md)

## Penalized Selection Index (PSI)

See [documentation](https://github.com/MarcooLopez/SFSI/blob/master/inst/doc/SFSI-extdoc-SSI.pdf) for the penalized selection index.

<div id="Outline" />

## Outline
  * [1. Data](#data)    
  * [2. Phenotypic vs Genotypic Selection Index](#GSI&PSI)
  * [3. Penalized Phenotypic and Genotypic Selection Index](#Penalized)
  * [4. Accuracy of the index](#SGSI&SPSIcv)
  * [5. Prediction of non-observed data](#predict)
  * [6. Genotypic covariance components patterns](#Scen)
   
-------------------------------------------------------------------------------------------

<div id="data" />

### 1. Data

Data will be simulated for *n* observations and *p* predictors. Both phenotypic values of the response variable *y* and predictors *X* are generated as the sum of some **breeding value** plus some environmental deviation in such a way that there is a given correlation between the phenotypic and breeding values (heritability). Also, some correlation exists between the breeding value of the response and that of all predictors (co-heritabilities), this value is equivalent to the squared root of the genetic correlation.

```r

# Function to simulate data
simulate_data <- function(n,p,h2y,h2xy,h2x,seed)
{
  set.seed(seed)
  
  # Response variable 
  Uy = sqrt(h2y)*as.vector(scale(rnorm(n)))    # Breeding value
  Ey =  sqrt(1-h2y)*as.vector(scale(rnorm(n))) # Environmental term
  y = scale(Uy + Ey)

  # Simulating predictor variables
  Ux <- Ex <- x <- matrix(NA,ncol=p,nrow=n)

  for(j in 1:p)
  {
    a1 = sqrt(h2xy[j])*scale(Uy)
    a2 =  sqrt(1-h2xy[j])*scale(rnorm(n))
    Ux[,j] <- sqrt(h2x[j])*scale(a1 + a2)    # Breeding value
    Ex[,j] <- sqrt(1-h2x[j])*scale(rnorm(n)) # Environmental term
    x[,j] <- scale(Ux[,j] + Ex[,j])
  }
  
  return(list(y=y,Uy=Uy,x=x,Ux=Ux))
}

# Simulate data with a low heritable response 
n <- 1500
p <- 500

h2y <- 0.25             # Heritability of the response
h2xy <- rbeta(p,10,10)  # Co-heritabilities of the response with predictors
h2x <- rbeta(p,10,10)   # Heritabilities of the predictors

dat <- simulate_data(n,p,h2y,h2xy,h2x,1234)
y <- dat$y
Uy <- dat$Uy
x <- dat$x
Ux <- dat$Ux
```

Genotypic covariances (between response and predictors) can be calculated from this simulated data by calculating the covariance between the simulated breeding values; however in a real situation, breeding values are **un-observable** and covariances must be estimated from variance components using linear mixed models either using replicates or multivariate models considering kinship relationship among individuals
```r
# Genotypic covariance between response and predictors 
gencov <- as.vector(cov(Ux,Uy))
 
# Phenotypic covariance between response and predictors 
phencov <- as.vector(cov(x,y))

plot(phencov,gencov)

# Phenotypic covariance matrix among predictors
Px <- var(x)
```
[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="GSI&PSI" />

### 2. Phenotypic vs Genotypic selection index (SI)

An index that uses phenotypic covariances will yield predictions that are the best in predicting phenotypic values; however this approach is not a good practice when the goal is selecting the best genotypes (judged by their breeding values). In the later case, using genotypic covariances seems to be more appropiate.

Code below will calculate SI that use both phenotypic (PSI) and genotypic covariances (GSI) 
```r
# Regression coefficients for the indices
B1 <- as.vector(solve(Px)%*%gencov)
B2 <- as.vector(solve(Px)%*%phencov)

# Fited values (selection indices)
GSI <- as.vector(x %*% B1)
PSI <- as.vector(x %*% B2)
```

When predicting the breeding values, we can see how good they were predicted using both PSI and GSI by evaluating by their Mean Squared Error of prediction in genotypic value prediction
```r
# MSE of the phenotypic SI
sum((Uy-PSI)^2)/n

# MSE of the phenotypic SI
sum((Uy-GSI)^2)/n

library(ggplot2)
rg <- range(c(Uy,GSI,PSI))
dat <- rbind(data.frame(Uy=Uy,yHat=PSI,SI="PSI"),data.frame(Uy=Uy,yHat=GSI,SI="GSI"))
ggplot(dat,aes(Uy,yHat,color=SI,group=SI)) + lims(x=rg,y=rg) + theme_bw() + 
   geom_abline(slope=1, intercept=0,linetype=2,size=.5,color="gray60") + 
   geom_point(shape=21,size=0.8) + labs(x="True BV",y="Predicted BV") + 
   theme(legend.justification=c(1,0),legend.position=c(0.99,0.01)) 
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/Figure1_SI.png" width="380">
</p>

[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="Penalized" />

### 3. Penalized Phenotypic and Genotypic selection index

The penalized index is obtained by imposing a penalization in the estimation of the regression coefficients. The penalization is given by the parameter `lambda`.

Code below will calculate penalized genotypic SI (PGSI) and penalized phenotypic SI (PPSI) for 100 values of the penalization parameter
```r
library(SFSI)

nLambda <- 100

# Genotypic SI
fm1 <- solveEN(Px,gencov,nLambda=nLambda)

# Phenotypic SI
fm2 <- solveEN(Px,phencov,nLambda=nLambda)

# Regression coefficients for each value of lambda
B1 <- as.matrix(fm1$beta)
B2 <- as.matrix(fm2$beta)

# Fitted values (selection indices)
PGSI <- x %*% t(B1)
PPSI <- x %*% t(B2)
```

The performance of the index is assessed by its accuracy of selection (i.e., correlation between the index and the breeding values) which  varies as the penalization parameter lambda changes

Again, the accuracy can be calculated using this simulated data but must be inferred from variance components in real data. Code below computes the accuracy of the index along the parameter lambda (logarithm scale)
```r
# Accuracy of the indices
accPGSI <- as.vector(cor(Uy,PGSI))
accPPSI <- as.vector(cor(Uy,PPSI))

dat <- rbind(
  data.frame(SI="PGSI",accuracy=accPGSI,df=fm1$df,lambda=fm1$lambda),
  data.frame(SI="PPSI",accuracy=accPPSI,df=fm2$df,lambda=fm2$lambda)
)

# Accuracy of the non-penalized SI (canonical SI)
dat2 <- data.frame(SI=c("PGSI","PPSI"),x=min(dat$lambda),accuracy=cor(cbind(Uy,GSI,PSI))[1,-1])

ggplot(dat[dat$df>1,],aes(-log(lambda),accuracy,color=SI,group=SI)) + 
    geom_line(size=0.8) + theme_bw() +
    geom_point(data=dat2,aes(-log(x),accuracy),size=2) 
```
[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="SGSI&SPSIcv" />

### 4. Accuracy of the penalized phenotypic vs penalized genotypic SI using cross-validation

The performance of the index is assessed by its accuracy of selection which varies according to the penalization parameter lambda.

Again, the accuracy can be calculated using this simulated data but must be inferred from variance components in real data. The accuracy is equal to the product of the squared root of the heritability of the index times the genetic correlation between the index and the target.

Code below will calculate the accuracy of the GSI and PSI in a cross validation fashion. The total number of data points will be `nRep` times `nFold`

```r
# Function to perform cross-validation  
SI_CV <- function(x,y,Ux,Uy,covType=c("geno","pheno"),nRep,nFold,nLambda,tol=1E-4,maxIter=500)
{
  n <- length(y)
  
  # Objects to store results
  accSI <- dfSI <- lambdaSI <- c()
  
  seeds <- seq(1E2,.Machine$integer.max,length=nRep)
  for(rep in 1:nRep)
  {
    # Create folds
    set.seed(seeds[rep])
    folds <- rep(seq(1:nFold), ceiling(n/nFold))[1:n]
    folds <- sample(folds)

    for(k in 1:nFold)
    {
      trn <- which(folds != k)
      tst <- which(folds == k)

      xTRN <- scale(x[trn,])
      
      # Covariances between response and predictors 
      if(covType=="pheno") covariance <- drop(cov(xTRN,scale(y[trn])))
      if(covType=="geno") covariance <- drop(cov(Ux[trn,],Uy[trn]))
 
      # Phenotypic covariance matrix among predictors
      Px <- var(xTRN)

      # Estimate regression coefficients
      fm <- solveEN(Px,covariance,tol=tol,maxIter=maxIter,nLambda=nLambda)
  
      # Retrieve data from 'df' and 'lambda'
      dfSI <- cbind(dfSI,fm$df)
      lambdaSI <- cbind(lambdaSI,fm$lambda)

      # Calculate the index (in testing set)
      SI <- fitted(fm,scale(x[tst,]))  
 
      # Accuracy of the index (correlation with the TARGET in testing set)
      corre <- suppressWarnings(cor(Uy[tst],SI))
      accSI <- cbind(accSI,as.vector(corre))
      cat("----- Rep=",rep,". Fold=",k,". Done\n")
    }
  }
  return(list(accSI=accSI,dfSI=dfSI,lambdaSI=lambdaSI))
}


nRep <- 10      # Number of replicates of CV
nFold <- 5     # Number of folds
nLambda <- 100   # Number of indices to generate

# Perform CV
out1 <- SI_CV(x,y,Ux,Uy,"geno",nRep,nFold,nLambda,tol=2E-4)     # GSI
out2 <- SI_CV(x,y,Ux,Uy,"pheno",nRep,nFold,nLambda,tol=2E-4)    # PSI

# Average accuracy of the indices across folds
accPGSI <- apply(out1$accSI,1,mean)
accPPSI <- apply(out2$accSI,1,mean)

dfPGSI <- apply(out1$dfSI,1,mean)
dfPPSI <- apply(out2$dfSI,1,mean)
lambdaPGSI <- apply(out1$lambdaSI,1,mean)
lambdaPPSI <- apply(out2$lambdaSI,1,mean)

dat <- rbind(
  data.frame(SI="PGSI",accuracy=accPGSI,df=dfPGSI,lambda=lambdaPGSI),
  data.frame(SI="PPSI",accuracy=accPPSI,df=dfPPSI,lambda=lambdaPPSI)
)

# Plot the average accuracy (in testing set) across all fold-replications
ggplot(dat[dat$df>1,],aes(-log(lambda),accuracy,color=SI,group=SI)) + 
   geom_line(size=0.8) + theme_bw() +
   theme(legend.justification=c(1,1),legend.position=c(0.99,0.99)) 
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/Figure2_SI.png" width="410">
</p>

An optimal index can be obtained such as the accuracy is maximum. Code below will take the index with maximum accuracy within each fold-replication. The penalized index is compared with the non-penalized (canonical) SI 
```r
dat <- rbind(
  data.frame(sp="Canonical",SI="GSI",accuracy=out1$accSI[nrow(out1$accSI),]),
  data.frame(sp="Penalized",SI="GSI",accuracy=apply(out1$accSI,2,max,na.rm=TRUE)),
  data.frame(sp="Canonical",SI="PSI",accuracy=out2$accSI[nrow(out2$accSI),]),
  data.frame(sp="Penalized",SI="PSI",accuracy=apply(out2$accSI,2,max,na.rm=TRUE))
)

dat2 <- aggregate(accuracy~SI+sp,dat,mean)
rg <- range(dat$accuracy)
ggplot(dat,aes(SI,accuracy,fill=SI)) + stat_boxplot(geom="errorbar",width=0.2) + 
  facet_wrap(~sp,scales="free_y") + geom_boxplot(width=0.5) + 
  theme_bw() + theme(legend.position="none")
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/Figure3_SI.png" height="325">
</p>

[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="predict" />

### 5. Prediction of non-observed data

The problem of prediction of data that have not been observed will be mimicked by predicting a newly simulated dataset using a separated training dataset from where a value of the penalization parameter will be chosen.

```r
# Data to train
dat <- simulate_data(1200,p,h2y,h2xy=rbeta(p,2,10),h2x=rbeta(p,2,10),123)

# New data
newdat <- simulate_data(300,p,h2y,h2xy=rbeta(p,2,10),h2x=rbeta(p,2,10),1234)

# Cross validation in the training data
out1 <- SI_CV(dat$x,dat$y,dat$Ux,dat$Uy,"geno",5,3,nLambda,tol=1E-4,maxIter=200)    

# Get the lambda that achieved maximum accuracy
lambda0 <- c()
for(j in 1:ncol(out1$lambda)) lambda0[j] <- out1$lambdaSI[which.max(out1$accSI[,j]),j]

# Use the lambda obtained to calculate an Penalized GSI to predict the new data
lambda <- mean(lambda0)

x <- scale(rbind(dat$x,newdat$x))   # Use all data from trn and new dataset
Px <- var(x)  
covariance <- drop(cov(dat$Ux,dat$Uy))  # Must be estimated using linear models
fm1 <- solveEN(Px,covariance,lambda=lambda)

yHat <- fitted(fm1,newdat$x)

cor(yHat,newdat$Uy)   # Accuracy

# Comparison with the canonical GSI (lambda=0)
fm2 <- solveEN(Px,covariance,lambda=0)
cor(fitted(fm2,newdat$x),newdat$Uy)   # Accuracy
```


[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="Scen" />

### 6. Different patterns of phenotypic and genotypic covariances

The only difference between GSI and PSI is that GSI uses genotypic covariances instead of phenotypic covariances between predictors and response. The phenotypic correlation depends on all the heritability of the response, the heritability of the predictors, and the genetic correlation between predictors and the response. Whenever the heritabilities are very high, both phenotypic and genotypic correlations are almost the same and the accuracy of the GSI and PSI are quite similar; however different extents of heritability and genetic correlations are likely to yield GSI and PSI with different accuracies. The example above was done considering moderately heritable predictors with moderate genetic correlation.

Code below will perform the same analysis for four different escenarios:
1. Both heritability and genetic correlation are high
2. Low heritability and high genotypic correlation
3. High heritability and low genotypic correlation
4. Both heritability and genetic correlation are low

```r
# Set up
n <- 1500
p <- 500
h2y <- 0.25            # Heritability of the response

# Co-heritabilities (squared root of the genetic correlation) (high and low)
h2xyH <- rbeta(p,30,1)  
h2xyL <- rbeta(p,2,10)

# Heritabilities of the predictors (high and low)
h2xH <- rbeta(p,30,1) 
h2xL <- rbeta(p,2,10) 

# Function for relevant plots
library(ggpubr)
makePlot <- function(h2xy,h2x,out1,out2)
{
  thm0 <- theme_bw() + theme(plot.title = element_text(size=10,hjust=0.5),
    axis.title.x=element_blank()) 
  
  # Histogram of h2 and genetic correlation
  dat <- rbind(data.frame(comp="h2",x=h2x),data.frame(comp="r",x=h2xy^2))
  pp1 <- ggplot(dat,aes(x,fill=comp)) + labs(title="Heritability and \n genetic correlation") +
     geom_histogram(alpha=0.5,aes(y = ..density..),position='identity') + thm0
  
  # Plot of accuracy
  dat <- rbind( data.frame(SI="PGSI",accuracy=apply(out1$accSI,2,max,na.rm=TRUE)),
                data.frame(SI="PPSI",accuracy=apply(out2$accSI,2,max,na.rm=TRUE))
  )

  dat2 <- aggregate(accuracy ~ SI, dat,mean)
  rg <- range(dat$accuracy)
  pp2 <- ggplot(dat,aes(SI,accuracy,fill=SI)) + stat_boxplot(geom="errorbar",width=0.2) + 
    geom_boxplot(width=0.5) + ylim(rg[1]-0.12*diff(rg),rg[2]) + labs(title="\nAccuracy of the index") +
    thm0 + scale_fill_manual(values=c("#E59F00","#56B4E9")) + theme(legend.position="none") +
    annotate("text",x=dat2$SI,y=rg[1]-0.1*diff(rg),label=sprintf("%.3f", dat2$accuracy)) 

  ggarrange(pp1,pp2,widths=c(0.54,0.46))
}

```
**Scenario 1**
```r
dat <- simulate_data(n,p,h2y,h2xyH,h2xH,1234)

# Perform CV
out1 <- SI_CV(dat$x,dat$y,dat$Ux,dat$Uy,"geno",10,5,100,tol=5E-4,maxIter=200)    # GSI
out2 <- SI_CV(dat$x,dat$y,dat$Ux,dat$Uy,"pheno",10,5,100,tol=5E-4,maxIter=200)   # PSI

makePlot(h2xyH,h2xH,out1,out2)
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/Figure4_SI.png" height="330">
</p>

**Scenario 2**
```r
dat <- simulate_data(n,p,h2y,h2xyH,h2xL,1234)

# Perform CV
out1 <- SI_CV(dat$x,dat$y,dat$Ux,dat$Uy,"geno",10,5,100,tol=5E-4,maxIter=200)    # GSI
out2 <- SI_CV(dat$x,dat$y,dat$Ux,dat$Uy,"pheno",10,5,100,tol=5E-4,maxIter=200)   # PSI

makePlot(h2xyH,h2xL,out1,out2)
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/Figure5_SI.png" height="330">
</p>


**Scenario 3**
```r
dat <- simulate_data(n,p,h2y,h2xyL,h2xH,1234)

# Perform CV
out1 <- SI_CV(dat$x,dat$y,dat$Ux,dat$Uy,"geno",10,5,100,tol=5E-4,maxIter=200)    # GSI
out2 <- SI_CV(dat$x,dat$y,dat$Ux,dat$Uy,"pheno",10,5,100,tol=5E-4,maxIter=200)   # PSI

makePlot(h2xyL,h2xH,out1,out2)
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/Figure6_SI.png" height="330">
</p>


**Scenario 4**
```r
dat <- simulate_data(n,p,h2y,h2xyL,h2xL,1234)

# Perform CV
out1 <- SI_CV(dat$x,dat$y,dat$Ux,dat$Uy,"geno",10,5,100,tol=5E-4,maxIter=200)    # GSI
out2 <- SI_CV(dat$x,dat$y,dat$Ux,dat$Uy,"pheno",10,5,100,tol=5E-4,maxIter=200)   # PSI

makePlot(h2xyL,h2xL,out1,out2)
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/Figure7_SI.png" height="330">
</p>
