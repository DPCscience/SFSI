% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSI.R
\name{plotPath}
\alias{plotPath}
\title{Coefficients path plot}
\usage{
plotPath(fm, Z = NULL, K = NULL, indexK = NULL, 
   tst = NULL, title = NULL, maxCor = 0.85)
}
\arguments{
\item{fm}{An object of the 'SSI' or 'SFI' class}

\item{Z}{Design matrix for the random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \ifelse{html}{\out{<b>G</b> = <b>K</b>}}{\eqn{\textbf{G = K}}{G = K}}; otherwise \ifelse{html}{\out{<b>G</b> = <b>Z K Z'</b>}}{\eqn{\textbf{G} = \textbf{Z K Z'}}{G = Z K Z'}} is used}

\item{K}{Kinship relationships matrix. This can be a name of a binary file where the matrix is stored. Only needed when a \code{fm} object of the class 'SFI' is provided}

\item{indexK}{Vector of integers indicating which columns and rows will be read when \code{K} is the name of a binary file. Default \code{indexK=NULL} will read the whole matrix}

\item{tst}{Vector of integers indicating which individuals are in testing set and must be contained in \code{fm$tst}. Default \code{tst=NULL} will consider the whole vector \code{fm$tst} to plot}

\item{title}{A 'character' type string for the plot title}

\item{maxCor}{Maximum correlation allowed for two different coefficients. A group of coeffcients with a correlation greater than \code{maxCor} are likely to overlap in the plot thus only one is kept}
}
\description{
Plot the coefficients evolution path along the regularization parameter
}
\examples{
  require(SFSI)
  data(wheat.E1) 
  data(wheat.G) # Genomic relationship matrix
  X <- apply(WL[[9]],2,function(x)tapply(x,Y[,"gid"],mean))
  y = tapply(Y[,"YLD"],Y[,"gid"],mean)
  index <- intersect(rownames(G),names(y))
  G = G[index,index]
  X = X[index,]
  y = as.vector(y[index])
  
  # Sparse phenotypic regression
  fm1 = lars2(var(X),cov(y,X))
  
  # Sparse family index
  fm2 = SFI(y,K=G,tst=1:20,trn=21:length(y))
  
  \dontrun{
  plotPath(fm1)
  plotPath(fm2,maxCor=0.6)
  plotPath(fm2,K=G,maxCor=0.6)
  }
  
  # Path plot for the first individual in testing set for the SFI
  plotPath(fm2,K=G,tst=fm2$tst[1])
}
\author{
Marco Lopez-Cruz (\email{lopezcru@msu.edu}) and Gustavo de los Campos
}
