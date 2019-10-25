
\name{SFI_CV}
\alias{SFI_CV}
\title{Sparse Family Index Cross-Validation}
\usage{
SFI_CV(G, y, h2 = 0.5, trn = 1:length(y), nFolds = 5,
  indexG = NULL, kernel = NULL, lambda = NULL,
  nLambda = 100, method = c("CD1", "CD2", "GBLUP"),
  alpha = 1, seed = 123, mc.cores = getOption("mc.cores", 2L),
  tol = 1e-05, maxIter = 800, name = NULL, verbose = TRUE)
}
\arguments{
\item{G}{Genetic relatedness matrix. This can be a name of a binary file where the matrix is storaged}

\item{y}{Response variable}

\item{h2}{Heritability of the response variable. Default is \code{h2=0.5}}

\item{trn}{Vector of integers indicating which individuals are to be consider for the cross-validation.
Default is \code{trn=1:length(y)} will consider all individuals}

\item{nFolds}{Number of non-overlaping folds in which the data is split into. Options 2,3,5,10}

\item{indexG}{Vector of integers indicating which columns and rows will be read when \code{G} is the name of a binary file.
Default \code{indexG=NULL} will read the whole matrix}

\item{kernel}{Kernel transformation to be applied to \code{G[trn,trn]}. See \code{help(kernel2)} for details. Default \code{kernel=NULL} (no kernel)}

\item{lambda}{Penalization parameter sequence vector. Default is \code{lambda=NULL}, in this case a decreasing grid of \code{'nLambda'} lambdas will be generated
starting from a maximum equal to
\ifelse{html}{\out{<center><font face="Courier">max(abs(G[trn,tst])/alpha)</font></center>}}{\deqn{\code{max(abs(G[trn,tst])/alpha)}}{max(abs(G[trn,tst])/alpha)}}
to a minimum equal to zero. If \code{alpha=0} the grid is generated starting from a maximum equal to 5}

\item{nLambda}{Number of lambdas generated when \code{lambda=NULL}}

\item{method}{One of:
\itemize{
 \item \code{'CD1'}: Computes the coefficients for a provided grid of lambdas common to all individuals in testing set.
 \item \code{'CD2'}: Similar to \code{'CD1'} but using a grid of lambdas specific to each individual in testing set.
 \item \code{'GBLUP'}: Coefficients are derived with no penalization and they correspond to those of the kinship-based BLUP
}}

\item{alpha}{Numeric between 0 and 1 indicating the weights given to the L1 and L2-penalties}

\item{mc.cores}{Number of cores used to run the analysis in parallel. Default is \code{mc.cores=2}}

\item{tol}{Maximum error between two consecutive solutions of the iterative algorithm to declare convergence}

\item{maxIter}{Maximum number of iterations to run at each lambda step before convergence is reached}

\item{seed}{Numeric seed to fix randomization when creating folds}

\item{name}{Name given to the output for tagging purposes. Default \code{name=NULL} will give the name of the method used}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing each step}
}
\value{
List object containing the elements:
\itemize{
  \item \code{folds}: matrix containing the folds used for the cross-validation.
  \item \code{correlation}: matrix with the correlation between observed and predicted values (in testing set) within each fold (in rows).
  \item \code{accuracy}: matrix with the accuracy (in testing set) within each fold (in rows).
  \item \code{MSE}: matrix with the mean squared error of prediction (in testing set) within each fold (in rows).
  \item \code{lambda}: matrix with the sequence of values of lambda used (averaged across individuals) within each fold (in rows).
  \item \code{df}: matrix with the degrees of freedom (averaged across individuals) within each fold (in rows).
  \item \code{kernel}: transformation applied to the elements of \code{G}.
}
Elements used as inputs: \code{y}, \code{h2}, \code{trn}, \code{method}, \code{name}, are also returned. The returned object
is of the class 'SFI_CV' for which methods \code{plot} and \code{summary} exist
}
\description{
Fit the Sparse Family Index in a cross validation fashion by randomly splitting data into 2,3,4,5, or 10 folds.
These partitions are equivalent to 50:50, 66:33, 75:25, 80:20, and 90:10 \% training-testing.
}
\details{
Data is divided into \eqn{k}{k} folds and the SFI is sequentially calculated for (all individuals in) one fold using the remaining folds as 'training' set.

Solutions for an specfic fold are found using the 'SFI' function (see \code{help(SFI)}) with parameter \code{tst} indicating the individuals in that fold and \code{trn} indicating the individuals in the remaining folds.
}
\examples{
  \dontrun{
  require(SFSI)
  # Read data from BGLR package
  data(wheat,package="BGLR")
  X = scale(wheat.X)
  G = tcrossprod(X)/ncol(X)    # Genomic relationship matrix
  y = wheat.Y[,1]              # Response variable

  # Training and testing sets
  set.seed(1234)
  n = length(y)
  pTST = 0.3      # percentage to predict
  tst = sample(1:n,ceiling(pTST*n))
  trn = (1:n)[-tst]

  # Calculate heritability
  fm = BGLR::BGLR(y,ETA=list(list(K=G,model="RKHS")),nIter=10000,burnIn=3000,verbose=FALSE)
  varU = fm$ETA[[1]]$varU
  varE = fm$varE
  h2 = varU/(varU + varE)

  # Obtain a value of lambda via cross validation in training set
  fm1 = SFI_CV(G,y,h2,trn)
  plot(fm1)
  lambda = summary(fm1)$opt$lambda
  fm2 = SFI(G,y,h2,trn,tst,lambda=lambda)
  summary(fm2)$correlation        # Testing set accuracy

  fm3 = SFI(G,y,h2,trn,tst,lambda=0)
  summary(fm3)$correlation        # Testing set accuracy of the non-sparse index
  }
}
\references{
\itemize{
\item \insertRef{Efron2004}{SFSI}
\item \insertRef{Friedman2007}{SFSI}
\item \insertRef{Lush1947}{SFSI}
\item \insertRef{Perez2014}{SFSI}
\item \insertRef{VanRaden2008}{SFSI}
}
}
\author{
Marco Lopez-Cruz (\email{lopezcru@msu.edu}) and Gustavo de los Campos
}
\keyword{SFI_CV}
