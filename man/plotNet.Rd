% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSI.R
\name{plotNet}
\alias{plotNet}
\title{Network plot}
\usage{
plotNet(fm, B, Z = NULL, K, indexK = NULL, subsetG = NULL,
    tst = NULL, U = NULL, d = NULL, group = NULL, group.fill = NULL,
    group.shape = NULL, group.size = NULL, df = NULL, title,
    curve = FALSE, TST.col = "yellow", bg.col = "gray20", 
    line.col = "gray90", line.tick = 0.3, legend.pos = "topright",
    set = c("Testing", "Active", "Non-active"))
}
\arguments{
\item{fm}{An object of the 'SFI' class}

\item{B}{Matrix of regression coefficients with number of rows and number of columns equal to the length of the vectors provided in \code{fm$tst} and \code{fm$trn}, respectively}

\item{Z}{Design matrix for the random effects. When \code{Z=NULL} an identity matrix is considered (default) thus \code{G = Z}; otherwise \code{G = Z\%*\%K\%*\%t(Z)} is used}

\item{K}{Kinship relationships matrix. This can be a name of a binary file where the matrix is storaged}

\item{indexK}{Vector of integers indicating which columns and rows will be read when \code{K} is the name of a binary file. Default \code{indexK=NULL} will read the whole matrix}

\item{subsetG}{Vector of integers indicating which columns (and rows) from \code{G} the vectors \code{fm$trn} and \code{fm$tst} refer to. Default \code{subsetG=NULL} considers that elements \code{fm$trn} and \code{fm$tst} refer to columns (and rows) from \code{G}; otherwise elements in training and testing in \code{G} have indices \code{subsetG[fm$trn]}  and \code{subsetG[fm$tst]}}

\item{tst}{Vector of integers indicating which individuals are in testing set and must be contained in \code{fm$tst}. Default \code{tst=NULL} will consider the whole vector \code{fm$tst} to plot}

\item{U}{Matrix of principal components or eigenvectors from a spectral decomposition of \code{G} (see \code{help(eigen)})}

\item{d}{Vector of eigenvalues corresponding to each eigenvector}

\item{group}{Object of the 'dataframe' class with one column grouping for the individuals. The rows must match with the rows in \code{G} matrix}

\item{df}{Average number of individuals in the training set contributing to the prediction (active) of individuals in the testing set. Default \code{df=NULL} will use the \code{df} that yielded the optimal accuracy}

\item{title}{A 'character' type string for the plot title}

\item{TST.col}{A 'character' type string indicating the color cycle of points corresponding to individuals in testing set. Used only when \code{group} is provided}

\item{bg.col}{A 'character' type string indicating the plot background color}

\item{line.col}{A 'character' type string indicating the color of lines connecting 'active' training individuals with each individual in testing set}

\item{line.tick}{A numeric value indicating the tick of lines connecting 'active' training individuals with each individual in testing set}

\item{group.fill}{Vector of 'character' strings indicating the fill colors of each level of the grouping column provided \code{group}}

\item{group.shape}{Vector of 'integer' numbers indicating the shape of each level of the grouping column provided \code{group}}

\item{group.size}{Vector of 'numeric' values indicating the size of 'testing', 'active', 'non-active' and elements, respectively}

\item{curve}{\code{TRUE} or \code{FALSE} to whether draw curve of rect lines connecting 'active' training individuals with each individual in testing set}

\item{legend.pos}{Either "topright","bottomleft","bottomright","topleft", or "none" indicating where the legend is positioned in the plot}

\item{set}{Vector of 'character' strings indicating the names of the sets: testing group, predictors with non-zero coefficient, and predictors with zero coefficient in the SFI, respectively}
}
\description{
Plot for a set of individuals in a testing set, a network of individuals in training set in a SFI
}
\examples{
  require(SFSI)
  data(wheat.E1) 
  data(wheat.G) # Genomic relationship matrix
  y = tapply(Y[,"YLD"],Y[,"gid"],mean)
  index <- intersect(rownames(G),names(y))
  G = G[index,index]
  y = as.vector(y[index])
  
  fm = SFI(y,K=G,tst=1:20,trn=21:length(y))
  
  # Basic setting
  plotNet(fm,K=G)
  plotNet(fm,K=G,bg.col="white",line.col="gray25")
  
  # Passing a matrix of coeeficients
  B=as.matrix(coef(fm,df=15))
  plotNet(fm,B=B,K=G,curve=TRUE,group.size=c(3.5,1.5,1))
  
  # Using Spectral Value Decomposition and grouping
  EVD <- eigen(G)
  gp <- data.frame(group=kmeans(EVD$vectors[,1:2],centers=5)$cluster)
  plotNet(fm,curve=TRUE,group=gp,U=EVD$vectors,d=EVD$values)
}
\author{
Marco Lopez-Cruz (\email{lopezcru@msu.edu}) and Gustavo de los Campos
}
