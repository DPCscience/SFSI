# PFSI
## Penalized Family and Selection Indices

Prediction of **breeding values** for a target trait (*y*<sub>i</sub>) is usually done using indirect information:
1. Correlated traits measured in the same candidates
2. Measurements on the same trait of interest collected on related individuals

The first case correspond to a **Selection Index** while the second yield a **Family Index**.
All the available observation contribute to the predicted breeding value ![](https://latex.codecogs.com/gif.latex?u_i) of the *i*<sup>th</sup> candidate of selection as a linear combination of the form:
![](https://latex.codecogs.com/gif.latex?u_i%3D%5Ctextbf%7Bx%7D%27%5Cboldsymbol%7B%5Cbeta%7D_i),
where ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bx%7D) can be either some correlated traits measured in the same candidate, ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7Bx%7D_i%3D%28x_%7Bi1%7D%2C...%2Cx_%7Bip%7D%29%27), or measurements on the same trait collected on related individuals, ![](https://latex.codecogs.com/gif.latex?%5Ctextbf%7By%7D%3D%28y_%7B1%7D%2C...%2Cy_%7Bn%7D%29%27). 

The weights ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cbeta%7D_i%3D%28%5Cbeta_%7Bi1%7D%2C...%2C%5Cbeta_%7Bip%7D%29%27) are derived by minimizing the optimization problem:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Chat%7B%5Cbeta%7D%7D_i%3Darg%20min%5Cfrac%7B1%7D%7B2%7DE%5Cleft%28u_i-%5Ctextbf%7Bx%7D%27%5Cboldsymbol%7B%5Cbeta%7D_i%20%5Cright%29%5E2">
</p>
