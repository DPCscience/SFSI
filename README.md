# SFSI: An R-package to estimate sparse selection indices

The SFSI Package implements shrinkage and variable selection regression procedures into the Selection Indices framework. In this repository we maintain the latest version beta version. This is an extended version of the SFSI that contains data used in [Lopez-Cruz et al. (2020)](https://www.biorxiv.org/content/10.1101/625251v2) for the development of penalized selection indices.

### Family and Selection Indices

Prediction of **breeding values** (<img src="https://render.githubusercontent.com/render/math?math=u_i">) for a target trait (<img src="https://render.githubusercontent.com/render/math?math=y_i">) is usually done using a **Selection Index**. The prediction is done through indirect information:
1. Correlated traits measured in the same candidates
2. Measurements on the same trait of interest collected on related individuals

The second case is also refered to as **Family Index** since the borrowing of information is taken from genetic relateness among individuals. In the selection index all the available observation contribute to the prediction of the *i*<sup>th</sup> candidate of selection as a linear combination of the form: 
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=I_i=x_{i1}\beta_{i1} %2B x_{i2}\beta_{i2} %2B ... %2B x_{ip}\beta_{ip}">
</p>
or (in matrix notation)
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=I_i=\textbf{x}_i^t\boldsymbol{\beta}_i">
</p>

where the predictors 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{x}_i">
can be either some correlated traits measured in the same candidate, 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{x}_i=(x_{i1},...,x_{ip})^t">
, or measurements on the same trait collected on related individuals, 
<img src="https://render.githubusercontent.com/render/math?math=\textbf{y}=(y_1,...,y_p)^t">
. 

The weights 
![](https://render.githubusercontent.com/render/math?math=\boldsymbol{\beta}_i=(\beta_{i1},...,\beta_{ip})^t)
are derived by minimizing the optimization problem:
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\hat{\boldsymbol{\beta}}_i=\text{arg min}\frac{1}{2}E\left(u_i-\textbf{x}_i^t\boldsymbol{\beta}_i\right)^2">
</p>

Under standard assumptions, the solution to the above problem is 
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\hat{\boldsymbol{\beta}}_i=\textbf{P}_x^{-1}\textbf{G}_{xy}">
</p>

where <img src="https://render.githubusercontent.com/render/math?math=\textbf{P}_x"> is the phenotypic variance-covariance matrix among predictors <img src="https://render.githubusercontent.com/render/math?math=\textbf{x}_i">,  and <img src="https://render.githubusercontent.com/render/math?math=\textbf{G}_{xy}"> is the genetic covariances between predictors <img src="https://render.githubusercontent.com/render/math?math=\textbf{x}_i"> and response <img src="https://render.githubusercontent.com/render/math?math=y_i">.

### Penalized Indices
The regression coefficients can be derived by impossing a penalization in the above optimization function as
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\hat{\boldsymbol{\beta}}_i=\text{arg min}\left[\frac{1}{2}E\left(u_i-\textbf{x}_i^t\boldsymbol{\beta}_i\right)^2 %2B \lambda J(\boldsymbol{\beta}_i)\right]">
</p>

where 
<img src="https://render.githubusercontent.com/render/math?math=\lambda">
is a penalty parameter (![](https://render.githubusercontent.com/render/math?math=\lambda=0)  yields the coefficients for the un-penalized index) and 
<img src="https://render.githubusercontent.com/render/math?math=J(\boldsymbol{\beta}_i)">
is a penalty function on the regression coefficients. Commonly used penalty functions are based on the L1 and L2 norms, 
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=L1:J(\boldsymbol{\beta}_i)=\sum_{j=1}^p{|\beta_{ij}|} \quad\quad L2:J(\boldsymbol{\beta}_i)=\frac{1}{2}\sum_{j=1}^p{\beta_{ij}^2}">
</p>

### Elastic-Net Penalized Index
An elastic-net penalized index considers a penalization being a weighted sum of both norms,
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=J(\boldsymbol{\beta}_i)=\alpha\sum_{j=1}^p{|\beta_{ij}|} %20%2B%20\frac{1}{2}(1-\alpha)\sum_{j=1}^p{\beta_{ij}^2}">
</p>

where <img src="https://render.githubusercontent.com/render/math?math=\alpha"> is a weighting parameter. Therefore the optimization problem becomes

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\hat{\boldsymbol{\beta}}_i=\text{arg min}\left[\frac{1}{2}E\left(u_i-\textbf{x}_i^t\boldsymbol{\beta}_i\right)^2 %2B \lambda \alpha\sum_{j=1}^p{|\beta_{ij}}| %20%2B%20\frac{1}{2}\lambda(1-\alpha)\sum_{j=1}^p{\beta_{ij}^2}\right]">
</p>

The L1-penalized and L2-penalized indices appear as special cases of the Elastic-Net-penalized index when
<img src="https://render.githubusercontent.com/render/math?math=\alpha=1">
 and
<img src="https://render.githubusercontent.com/render/math?math=\alpha=0">
 , respectively. When <img src="https://render.githubusercontent.com/render/math?math=\alpha=0">, the solution has closed form:

<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\hat{\boldsymbol{\beta}}_i=\left(\textbf{P}_x%2B\lambda\textbf{I}\right)^{-1}\textbf{G}_{xy}">
</p>

If <img src="https://render.githubusercontent.com/render/math?math=\alpha \gt 0">
, no closed form solution exists; however, a solution can be obtained using iterative algorithms such as Least Angle Regression (LARS) (Efron, 2004) or Coordinate Descent algorithms (Friedman, 2007).

### Sparse Family and Selection Indices using the SFSI R-package
Depending of the type of information used as predictors (either correlated traits measured in the same candidates or measurements on the same trait collected on related individuals), the problem can be seen either as a **Selection Index** or a **Family Index**. 
The penalized indices can be solved using the package SFSI that implements LARS and Coordinate Descent algorithms using as inputs <img src="https://render.githubusercontent.com/render/math?math=\textbf{P}_x"> and <img src="https://render.githubusercontent.com/render/math?math=\textbf{G}_{xy}">. The coefficients of the index are calculated for different values of <img src="https://render.githubusercontent.com/render/math?math=\lambda"> for a given value of the parameter
<img src="https://render.githubusercontent.com/render/math?math=\alpha">
. Optimal indices can be obtained by choosing the values of these parameters that maximize the accuracy.

### Documentation
* **[Penalized Selection Index](https://github.com/MarcooLopez/PFSI/blob/master/inst/md/selection_index.md)**
* **[Sparse Family Index](https://github.com/MarcooLopez/PFSI/blob/master/inst/md/family_index.md)**

**Package installation from Github**

Installation of SFSI package requires a R-version greater than 3.5.0
```r
  install.packages('devtools',repos='https://cran.r-project.org/')      #1. install devtools
  library(devtools)                                                     #2. load the library
  install_git('https://github.com/MarcooLopez/SFSI')                    #3. install SFSI from GitHub
```

**How to cite SFSI R-package:**
* Lopez-Cruz M, Olson E, Rovere G, Crossa J, Dreisigacker S, Mondal S, Sing R & de los Campos G **(2020)**. Regularized selection indices for breeding value prediction using hyper-spectral image data. *Scientific Reports (in press)*.


### References
* Efron B, Hastie T, Johnstone I & Tibshirani R **(2004)**. Least angle regression. *The Annals of Statistics*, 32(2), 407–499.
* Friedman J, Hastie T, Höfling H & Tibshirani R **(2007)**. Pathwise coordinate optimization. *The Annals of Applied Statistics*, 1(2), 302–332.
