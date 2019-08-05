# PFSI
## Family and Selection Indices

Prediction of **breeding values** for a target trait (*y*<sub>*i*</sub>) is usually done using indirect information:
1. Correlated traits measured in the same candidates
2. Measurements on the same trait of interest collected on related individuals

The first case corresponds to what is called a **Selection Index** while the second yield a **Family Index**.
All the available observation contribute to the predicted breeding value *u*<sub>*i*</sub> of the *i*<sup>th</sup> candidate of selection as a linear combination of the form:
![](https://latex.codecogs.com/gif.latex?u_i%3D%5Ctextbf%7Bx%7D%27%5Cboldsymbol%7B%5Cbeta%7D_i),
where the predictors ***x*** can be either some correlated traits measured in the same candidate, ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bx%7D_i%3D%28x_%7Bi1%7D%2C...%2Cx_%7Bip%7D%29), or measurements on the same trait collected on related individuals, ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D%3D%28y_1%2C...%2Cy_p%29). 

The weights ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cbeta%7D_i%3D%28%5Cbeta_%7Bi1%7D%2C...%2C%5Cbeta_%7Bip%7D%29) are derived by minimizing the optimization problem:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Chat%7B%5Cbeta%7D%7D_i%3D%5Ctext%7Barg%20min%7D%5Cfrac%7B1%7D%7B2%7DE%5Cleft%28u_i-%5Ctextbf%7Bx%7D%27%5Cboldsymbol%7B%5Cbeta%7D_i%5Cright%29%5E2">
</p>

Under standard assumptions, the solution to the above problem is 
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Chat%7B%5Cbeta%7D%7D_i%3D%5Ctextbf%7BP%7D%5E%7B-1%7D_x%5Ctextbf%7BG%7D_%7Bxy%7D">
</p>

where ***P***<sub>*x*</sub> is the phenotypic variance-covariance matrix among ***x*** and ***G***<sub>*xy*</sub> is the genetic covariances between predictors ***x*** and response *y*.

## Penalized Indices
The regression coefficients can be derived by impossing a penalization in the above optimization function as
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Chat%7B%5Cbeta%7D%7D_i%3D%5Ctext%7Barg%20min%7D%5Cleft%5B%5Cfrac%7B1%7D%7B2%7DE%5Cleft%28u_i-%5Ctextbf%7Bx%7D%27%5Cboldsymbol%7B%5Cbeta%7D_i%5Cright%29%5E2&plus;%5Clambda%20J%28%5Cboldsymbol%7B%5Cbeta%7D_i%29%5Cright%5D">
</p>

where ![](https://latex.codecogs.com/gif.latex?%5Clambda) is a penalty parameter (![](https://latex.codecogs.com/gif.latex?%5Clambda%3D0) yields the coefficients for the un-penalized index) and ![](https://latex.codecogs.com/gif.latex?J%28%5Cboldsymbol%7B%5Cbeta%7D%29) is a penalty function. Commonly used penalty functions are based on the L1 and L2 norms, 
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?L1%3AJ%28%5Cboldsymbol%7B%5Cbeta%7D%29%3D%5Csum_%7Bj%3D1%7D%5Ep%7C%5Cbeta_j%7C%5Cqquad%20%5Cqquad%20L2%3A%20J%28%5Cboldsymbol%7B%5Cbeta%7D%29%3D%5Cfrac%7B%7D%7B1%7D%7B2%7D%5Csum_%7Bj%3D1%7D%5Ep%5Cbeta_j%5E2">
</p>

### Elastic-Net Penalized Index
An elastic-net penalized index considers a penalization being a weighted sum of both norms,
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?J%28%5Cboldsymbol%7B%5Cbeta%7D%29%3D%5Calpha%5Csum_%7Bj%3D1%7D%5Ep%7C%5Cbeta_j%7C&plus;%5Cfrac%7B1%7D%7B2%7D%281-%5Calpha%29%5Csum_%7Bj%3D1%7D%5Ep%5Cbeta_j%5E2">
</p>

Therefore the optimization problem becomes,
<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Chat%7B%5Cbeta%7D%7D_i%3D%5Ctext%7Barg%20min%7D%5Cleft%5B%5Cfrac%7B1%7D%7B2%7DE%5Cleft%28u_i-%5Ctextbf%7Bx%7D%27%5Cboldsymbol%7B%5Cbeta%7D_i%5Cright%29%5E2&plus;%5Clambda%20%5Calpha%5Csum_%7Bj%3D1%7D%5Ep%7C%5Cbeta_%7Bij%7D%7C&plus;%5Cfrac%7B1%7D%7B2%7D%5Clambda%281-%5Calpha%29%5Csum_%7Bj%3D1%7D%5Ep%5Cbeta_%7Bij%7D%5E2%29%5Cright%5D">
</p>

The L1-penalized and L2-penalized indices appear as special cases of the Elastic-Net-penalized index when ![](https://latex.codecogs.com/gif.latex?%5Calpha%3D1) and ![](https://latex.codecogs.com/gif.latex?%5Calpha%3D0), respectively. When ![](https://latex.codecogs.com/gif.latex?%5Calpha%3D0), the solution has closed form:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Chat%7B%5Cboldsymbol%7B%5Cbeta%7D%7D_i%3D%5Cleft%28%5Ctextbf%7BP%7D_x&plus;%5Clambda%5Ctextbf%7BI%7D%20%5Cright%20%29%5E%7B-1%7D%5Ctextbf%7BG%7D_%7Bxy%7D">
</p>

If ![](https://latex.codecogs.com/gif.latex?%5Calpha%3E0), no closed form solution exists; however, a solution can be obtained using iterative algorithms such as Least Angle Regression (LARS) (Efron, 2004) or Coordinate Descent algorithms (Friedman, 2007).

### Penalized Family and Selection Indices using 'PFSI' R-package
Depending of the type of information used as predictors ***x*** (either correlated traits measured in the same candidates or measurements on the same trait collected on related individuals), the problem can be seen either as a Selection or a Family Index. 
The penalized indices can be solved using the package 'PFSI' that implements LARS and Coordinate Descent algorithms using as inputs ***P***<sub>*x*</sub> and ***G***<sub>*xy*</sub>. The coefficients of the index are calculated for different values of ![](https://latex.codecogs.com/gif.latex?%5Clambda) for a given value of the parameter ![](https://latex.codecogs.com/gif.latex?%5Calpha). Optimal indices can be obtained by choosing the values of these parameters that maximize the accuracy.

**Package installation from Github**
```r
  install.packages('devtools',repos='https://cran.r-project.org/')      #1. install devtools
  library(devtools)                                                     #2. load the library
  install_git('https://github.com/MarcooLopez/PFSI')                    #3. install PFSI from GitHub
```

* **[Penalized Selection Index](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/single_environment.md)**
* **[Penalized Family Index](https://github.com/MarcooLopez/Genomic-Selection-Demo/blob/master/multi_environment.md)**
