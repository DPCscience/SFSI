
### 1. Download and prepare data
```r
library(BGLR)
data(wheat)
A <- wheat.A
X <- wheat.X
Y <- wheat.Y

# Calculating the genomic relationship matrix
G <- tcrossprod(scale(X))/ncol(X)
```
