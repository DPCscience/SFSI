* [Back to main page](https://github.com/MarcooLopez/SFSI/blob/master/README.md)

## Sparse Selection Index (SSI)

<div id="Outline" />

## Outline
  * [1. Data](#data)    
  * [2. Phenotypic vs Genotypic Selection Index](#GSI&PSI)
  * [3. Accuracy of the index](#GSI&PSIcv)
  * [4. Genotypic covariance components patterns](#Scen)
   
-------------------------------------------------------------------------------------------

<div id="data" />

### 1. Data

Data will be simulated for *n* observations and *p* predictors. Both phenotypic values of the response variable *y* and predictors *X* are generated as the sum of some genotypic value plus some environmental deviation in such a way that there is a given correlation between the phenotypic and genotypic values (heritability). Also, some correlation exists between genotypic value of the response and that of all predictors (co-heritabilities), this value is equivalent to the squared root of the genetic correlation.

```r

# Function to simulate data
simulate_data <- function(n,p,h2y,h2xy,h2x,seed)
{
  set.seed(seed)
  
  # Response variable 
  Uy = sqrt(h2y)*as.vector(scale(rnorm(n)))    # Genotypic value
  Ey =  sqrt(1-h2y)*as.vector(scale(rnorm(n))) # Environmental term
  y = scale(Uy + Ey)

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
  
  return(list(y=y,Uy=Uy,x=x,Ux=Ux))
}

# Simulate data with a low heritable response 
n <- 1500
p <- 750

h2y <- 0.25          # Heritability of the response
h2xy <- rbeta(p,10,2)  # Co-heritabilities of the response with predictors
h2x <- rbeta(p,10,10)   # Heritabilities of the predictors

dat <- simulate_data(n,p,h2y,h2xy,h2x,1234)
y <- dat$y
Uy <- dat$Uy
x <- dat$x
Ux <- dat$Ux
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
[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="GSI&PSI" />

### 2. Phenotypic vs Genotypic selection index

An index that uses phenotypic covariances will yield predictions that are the best in predicting phenotypic values; however this approach is not a good practice when the goal is selecting the best genotypes (judged by their genotypic values). In the later case, using genotypic covariances seems to be more appropiate.

Code below will calculate SS that use both phenotypic (PSI) and genotypic covariances (GSI) 
```r
library(SFSI)

# Genotypic SS
fm1 <- SSI(Px,gencov,method="CD")

# Phenotypic SS
fm2 <- SSI(Px,phencov,method="CD")

# Regression coefficients for each value of lambda
B1 <- as.matrix(fm1$beta)
B2 <- as.matrix(fm2$beta)

# Fited values (selection indices)
GSI <- (x %*% t(B1))
PSI <- (x %*% t(B2))
```

The resulting indices are obtained such that the correlation between the index and the target is maximum (accuracy of selection). In this cases, the target of the phenotypic SS is the phenotype and for the genotypic SS is the genotype.
The indices can be also evaluated by their Mean Squared Error of prediction in genotypic value prediction
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
[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="GSI&PSIcv" />

### 3. Phenotypic vs Genotypic selection index using cross-validation

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
  
  for(rep in 1:nRep)
  {
    # Create folds
    set.seed(rep*1234)
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
      fm <- SSI(Px,covariance,method="CD",tol=tol,maxIter=maxIter,nLambda=nLambda)
  
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
out1 <- SI_CV(x,y,Ux,Uy,"geno",nRep,nFold,nLambda)   # GSI
out2 <- SI_CV(x,y,Ux,Uy,"pheno",nRep,nFold,nLambda)    # PSI

# Accuracy of the indices across folds
accGSI <- apply(out1$accSI,1,mean)
accPSI <- apply(out2$accSI,1,mean)

dfGSI <- apply(out1$dfSI,1,mean)
dfPSI <- apply(out2$dfSI,1,mean)
lambdaGSI <- apply(out1$lambdaSI,1,mean)
lambdaPSI <- apply(out2$lambdaSI,1,mean)

dat <- rbind(
  data.frame(SI="GSI",accuracy=accGSI,df=dfGSI,lambda=lambdaGSI),
  data.frame(SI="PSI",accuracy=accPSI,df=dfPSI,lambda=lambdaPSI)
)

# Plot the average accuracy (in testing set) across all fold-replications
ggplot(dat[dat$df>1,],aes(-log(lambda),accuracy,color=SI,group=SI)) + 
   geom_line(size=0.8)

```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/CV_lambda_1.png" width="415">
</p>

An optimal index can be obtained such as the accuracy is maximum. Code below will take the index with maximum accuracy within each fold-replication 
```r
dat <- rbind(
  data.frame(SI="GSI",accuracy=apply(out1$accSI,2,max,na.rm=TRUE)),
  data.frame(SI="PSI",accuracy=apply(out2$accSI,2,max,na.rm=TRUE))
)

rg <- range(dat$accuracy)
ggplot(dat,aes(SI,accuracy,fill=SI)) + stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(width=0.5) + ylim(rg[1]*0.99,ifelse(rg[2]*1.01>1,1,rg[2]*1.01))
  
```

<p align="center">
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/CV_mean_1.png" width="415">
</p>

[Back to Outline](#Outline)

-------------------------------------------------------------------------------------------

<div id="Scen" />

### 4. Different patterns of phenotypic and genotypic covariances

The only difference between GSI and PSI is that GSI uses genotypic covariances instead of phenotypic covariances between predictors and response. Thus, whenever they are very different the results of GSI and PSI are expected to be different. The genotypic covariances will depend of the heritability of the predictors and the genetic correlation between predictors and the response. The example above was done considering moderately heritable predictors with high genetic correlation.

Code below will perform the same analysis for four different escenarios:
1. Both heritability and genetic correlation are high
2. Low heritability and high genotypic correlation
3. High heritability and low genotypic correlation
4. Both heritability and genetic correlation are low

```r
# Set up
n <- 1500
p <- 750
h2y <- 0.25            # Heritability of the response

# Co-heritabilities (squared root of the genetic correlation) (high and low)
h2xyH <- rbeta(p,30,1)  
h2xyL <- rbeta(p,2,10)

# Heritabilities of the predictors (high and low)
h2xH <- rbeta(p,30,1) 
h2xL <- rbeta(p,2,10) 

# Funtion for relevant plots
library(ggpubr)
makePlot <- function(h2xy,h2x,out1,out2)
{
  thm0 <- theme(plot.title = element_text(hjust = 0.5))
  
  # Histogram of h2 and genetic correlation
  dat <- rbind(data.frame(comp="h2",x=h2x),data.frame(comp="r",x=h2xy^2))
  pp1 <- ggplot(dat, aes(x, fill = comp)) + labs(title="Heritability and \n genetic correlation") +
     geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') + thm0
  
  # Plot of accuracy
  dat <- rbind( data.frame(SI="GSI",accuracy=apply(out1$accSI,2,max,na.rm=TRUE)),
                data.frame(SI="PSI",accuracy=apply(out2$accSI,2,max,na.rm=TRUE))
  )

  dat2 <- aggregate(accuracy ~ SI, dat,mean)
  rg <- range(dat$accuracy)
  pp2 <- ggplot(dat,aes(SI,accuracy,fill=SI)) + stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(width=0.5) + ylim(rg[1]*0.995,ifelse(rg[2]*1.005>1,1,rg[2]*1.005)) +
    labs(title="Accuracy of the GSI and PSI")+theme(legend.position="none")+thm0+
    annotate("text", x=dat2$SI, y=rg[1]*0.996, label=sprintf("%.3f", dat2$accuracy))

  ggarrange(pp1,pp2)
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
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/CV_sce_1.png" height="355">
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
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/CV_sce_2.png" height="355">
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
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/CV_sce_3.png" height="355">
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
<img src="https://github.com/MarcooLopez/SFSI/blob/master/inst/md/CV_sce_4.png" height="355">
</p>
