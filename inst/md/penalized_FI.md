## Penalized Family Index
The kinship-based BLUPcanonical Family Index (i.e., un-penalized),

### 1. Data
Data from CIMMYT’s Global Wheat Program. Lines were evaluated for grain yield (each entry corresponds to an average of two plot records) at four different environments; phenotypes (*wheat.Y* object) were centered and standardized to a unit variance within environment. Each of the lines were genotyped for 1279 diversity array technology (DArT) markers. At each marker two homozygous genotypes were possible and these were coded as 0/1. Marker genotypes are given in the object *wheat.X*. Finally a matrix *wheat.A* provides the pedigree relationships between lines computed from the pedigree records. Data is available for download in the R-package 'BGLR'.

**Download data**
```r
library(BGLR)
data(wheat)
A <- wheat.A
X <- wheat.X
Y <- wheat.Y

# Calculating the genomic relationship matrix
G <- tcrossprod(scale(X))/ncol(X)
```

### 2. Comparing predictive ability with the un-penalized family index
