
\name{SFI}
\alias{SFI}
\title{Sparse Family Index}
\usage{
SFI(G, y, h2 = 0.5, trn = 1:length(y), tst = 1:length(y),
  indexG = NULL, subset = NULL, kernel = NULL, lambda = NULL,
  nLambda = 100, method = c("CD1", "CD2"), alpha = 1, name = NULL,
  mc.cores=getOption("mc.cores", 2L), tol = 1E-5, maxIter = 800,
  saveAt = NULL, verbose = TRUE)
}
\arguments{
\item{G}{Genetic relatedness matrix. This can be a name of a binary file where the matrix is storaged}

\item{y}{Response variable}

\item{h2}{Heritability of the response variable. Default is \code{h2=0.5}}

\item{trn}{Vector of integers indicating which individuals are in training set. Default is \code{trn=1:length(y)} will consider all individuals as training}

\item{tst}{Vector of integers indicating which individuals are in testing set. Default is \code{tst=1:length(y)} will consider all individuals as testing}

\item{indexG}{Vector of integers indicating which columns and rows will be read when \code{G} is the name of a binary file.
Default \code{indexG=NULL} will read the whole matrix}

\item{subset}{A two-elements numeric vector \eqn{c(m,M)} to fit the model only for the \ifelse{html}{\out{m<sup>th</sup>}}{\eqn{m^{th}}{m^{th}}}
subset out of \eqn{M} subsets that the
testing set will be divided into. Results can be automatically saved when \code{saveAt} parameter is provided and can be retrieved later
using function \code{collect}. Default is \code{subset=NULL} for no subsetting, then the model is fitted using all information}

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
}
Default is \code{method='CD1'}}

\item{name}{Name given to the output for tagging purposes. Default \code{name=NULL} will give the name of the method used}

\item{mc.cores}{Number of cores used to run the analysis in parallel. Default is \code{mc.cores=2}}

\item{tol}{Maximum error between two consecutive solutions of the iterative algorithm to declare convergence}

\item{maxIter}{Maximum number of iterations to run at each lambda step before convergence is reached}

\item{saveAt}{Prefix name that will be added to the output files name to be saved, this may include a path. Regression coefficients
are saved as binary file as a single-precision (32 bits, 7 significant digits) variable. Default \code{saveAt=NULL} will no save any output}

\item{verbose}{\code{TRUE} or \code{FALSE} to whether printing each step}

\item{alpha}{Numeric between 0 and 1 indicating the weights given to the L1 and L2-penalties}
}
\value{
List object containing the elements:
\itemize{
  \item \code{BETA}: list object containing, for each individual in testing set, a matrix of regression coefficients.
  \item \code{alpha}: value for the elastic-net weights used.
  \item \code{lambda}: matrix with the sequence of values of lambda used (for each individual in rows).
  \item \code{df}: degrees of freedom (averaged across individuals), number of non-zero predictors at each solution.
  \item \code{kernel}: transformation applied to the elements of \code{G}.
}
Elements used as inputs: \code{y}, \code{h2}, \code{trn}, \code{tst}, \code{method}, \code{name}, are also returned.
The returned object is of the class 'SFI' for which methods \code{fitted}, \code{plot} and \code{summary} exist
}
\description{
Computes the entire Elastic-Net solution for the regression coefficients of a Family Index simultaneously for all values of the penalization parameter via the Coordinate Descent (CD) algorithm.
}
\details{
Finds solutions for the genetic values of individuals in a linear model

\ifelse{html}{\out{<center>y<sub>i</sub> = u<sub>i</sub> + e<sub>i</sub></center>}}{\deqn{y_i=u_i+e_i}{y_i = u_i + e_i}}

where
\ifelse{html}{\out{y<sub>i</sub>}}{\eqn{y_i}{y_i}} is the response for the \ifelse{html}{\out{i<sup>th</sup>}}{\eqn{i^{th}}{i^th}} observation,
\ifelse{html}{\out{u<sub>i</sub>}}{\eqn{u_i}{u_i}}
is the genetic value of the individual, and
\ifelse{html}{\out{e<sub>i</sub>}}{\eqn{e_i}{e_i}}
is an environmental residual.

The values \ifelse{html}{\out{u<sub>i</sub>}}{\eqn{u_i}{u_i}} for
\ifelse{html}{\out{i = 1,2,...,n<sub>tst</sub>}}{\eqn{i=1,2,...,n_{tst}}{i = 1,2,...,n_tst}} individuals in a testing set are estimated using (as predictors) all available observations in a training set as

\ifelse{html}{\out{<center>u<sub>i</sub> = <b>y</b>'<sub>trn</sub> <b>&beta;</b><sub>i</sub></center>}}{\deqn{u_i=\textbf{y}_{trn}'\boldsymbol{\beta}_i}{u_i = y'_trn beta_i}}

where \ifelse{html}{\out{<b>&beta;</b><sub>i</sub>}}{\eqn{\boldsymbol{\beta}_i}{beta_i}}
is a vector of regression coefficients which are estimated as function of the genetic covariance between predictors
and that individual (\ifelse{html}{\out{<b>G</b><sub>trn,tst(i)</sub>}}{\eqn{\textbf{G}_{trn,tst(i)}}{G[trn,tst(i)]}}) and the genetic variance among predictors (\ifelse{html}{\out{<b>G</b><sub>trn,trn</sub>}}{\eqn{\textbf{G}_{trn,trn}}{G[trn,trn]}}).

These coefficients are obtained, separately for each individual in the testing set, by minimizing the penalized mean squared error function

\ifelse{html}{\out{<center>-<b>G</b>'<sub>trn,tst(i)</sub> <b>&beta;</b><sub>i</sub> + 1/2 <b>&beta;</b>'<sub>i</sub>(<b>G</b><sub>trn,trn</sub> + &lambda;<sub>0</sub><b>I</b>)<b>&beta;</b><sub>i</sub> + &lambda; J(<b>&beta;</b><sub>i</sub>)</center>}}{\deqn{-\textbf{G}_{trn,tst(i)}' \boldsymbol{\beta}_i + 1/2 \boldsymbol{\beta}_i'(\textbf{G}_{trn,trn} + \lambda_0\textbf{I}) \boldsymbol{\beta}_i + \lambda J(\boldsymbol{\beta}_i)}{-G'[trn,tst(i)] beta_i + 1/2 beta_i'(G[trn,trn] + lambda_0 I)beta_i + lambda J(beta_i)}}

where \ifelse{html}{\out{&lambda;<sub>0</sub> = (1-h<sup>2</sup>)/h<sup>2</sup>}}{\eqn{\lambda_0=(1-h^2)/h^2}{lambda_0 = (1-h^2)/h^2}} is a fixed parameter expressed in terms of the heritability, \ifelse{html}{\out{&lambda;}}{\eqn{\lambda}{lambda}} is the penalization parameter, and
\ifelse{html}{\out{J(<b>&beta;</b><sub>i</sub>)}}{\eqn{J(\boldsymbol{\beta}_i)}{J(beta_i)}}
is a penalty function given by
\ifelse{html}{\out{<center>1/2(1-&alpha;)||<b>&beta;</b><sub>i</sub>||<sub>2</sub><sup>2</sup> + &alpha;||<b>&beta;</b><sub>i</sub>||<sub>1</sub></center>}}{\deqn{1/2(1-\alpha)||\boldsymbol{\beta}_i||_2^2 + \alpha||\boldsymbol{\beta}_i||_1}{1/2(1-alpha)||beta_i||_2^2 + alpha||beta_i||_1}}

where \ifelse{html}{\out{0 &le; &alpha; &le; 1}}{\eqn{0\leq\alpha\leq 1}{0 <= alpha <= 1}}, and
\ifelse{html}{\out{||<b>&beta;</b><sub>i</sub>||<sub>1</sub> = &sum;|&beta;<sub>ij</sub>|}}{\eqn{||\boldsymbol{\beta}_i||_1=\sum|\beta_{ij}|}{||beta_i||_1 = sum(|beta_ij|)}} and
\ifelse{html}{\out{||<b>&beta;</b><sub>i</sub>||<sub>2</sub><sup>2</sup> = &sum;&beta;<sub>ij</sub><sup>2</sup>}}{\eqn{||\boldsymbol{\beta}_i||_2^2=\sum\beta_{ij}^2}{||beta_i||_2^2 = sum(beta_ij^2)}} are the L1 and (squared) L2-norms, respectively.

Each individual solution is found using the 'solveEN' function (see \code{help(solveEN)}) by setting the parameter \code{XtX} equal to
\ifelse{html}{\out{<b>G</b><sub>trn,trn</sub> + &lambda;<sub>0</sub><b>I</b>}}{\eqn{\textbf{G}_{trn,trn}+\lambda_0\textbf{I}}{G[trn,trn] + lambda_0 I}}
and \code{Xty} equal to
\ifelse{html}{\out{<b>G</b><sub>trn,tst(i)</sub>}}{\eqn{\textbf{G}_{trn,tst(i)}}{G[trn,tst(i)]}}
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

  # Sparse family index
  fm = SFI(G,y,h2,trn,tst)

  yHat = fitted(fm)      # Predicted values for each SFI
  corTST = cor(y[tst],yHat)[1,]  # Testing set accuracy or
  tmp = summary(fm)
  corTST = tmp$correlation
  tmp$opt                # Optimal accuracy
  dat = data.frame(df=tmp$df,logLambda=-log(tmp$lambda),corTST)
  plot(corTST ~ logLambda, data=dat[dat$df>1,]) # Accuracy along the regularization parameter
  plot(fm)                 # Same plot as above
  plot(fm,py="MSE")        # MSE vs regularization
  plot(fm,G=G)             # Individual's training subsets for the optimal index
  plot(fm,G=G,df=10)       # Individual's training subsets
  plot(fm,path=TRUE)       # Coefficients path plot
  plot(fm,path=TRUE,G=G)   # Coefficients path plot showing kinship
  }
}
\references{
\itemize{
\item \insertRef{Efron2004}{SFSI}
\item \insertRef{Friedman2007}{SFSI}
\item \insertRef{Hoerl1970}{SFSI}
\item \insertRef{Lush1947}{SFSI}
\item \insertRef{Perez2014}{SFSI}
\item \insertRef{Tibshirani1996}{SFSI}
\item \insertRef{VanRaden2008}{SFSI}
\item \insertRef{Zou2005}{SFSI}
}
}
\author{
Marco Lopez-Cruz (\email{lopezcru@msu.edu}) and Gustavo de los Campos
}
\keyword{SFI}
