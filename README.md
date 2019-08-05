# PFSI
## Penalized Family and Selection Indices

Prediction of **breeding values** for a target trait is usually done using indirect information:
1. Correlated traits measured in the same candidates
2. Measurements on the trait of interest collected on related individuals

The first case correspond to a **Selection Index** while the second yield a **Family Index**.
All the available observations (***x***) contribute to the predicted breeding value of each candidate of selection in a linear combination of the form:
![](https://latex.codecogs.com/gif.latex?u_i%3D%5Csum_%7Bk%3D1%7D%5Ep%7Bx_i%5Cbeta_%7Bik%7D%7D).

The weights ![](https://latex.codecogs.com/gif.latex?%5Cboldsymbol%7B%5Cbeta%7D_i%3D%28%5Cbeta_%7Bi1%7D%2C...%2C%5Cbeta_%7Bip%7D%29%27) are derived by minimizing the optimization problem

