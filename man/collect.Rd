% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_utils_wC.R
\name{collect}
\alias{collect}
\title{collect function}
\usage{
collect(prefix = "")
}
\arguments{
\item{prefix}{Prefix that was added to the output files name, this may include a path}
}
\value{
An object of the class 'SFI' for which methods \code{fitted}, \code{plot} and \code{summary} exist
}
\description{
Collects all outputs saved at the provided \code{saveAt} parameter from the SFI analysis when testing data was splited
according to \code{subset} parameter.
}
\examples{
  require(SFSI)
  data(wheat.E1) 
  data(wheat.G) # Genomic relationship matrix
  y = tapply(Y[,"YLD"],Y[,"gid"],mean)  # Response variable
  index <- intersect(rownames(G),names(y))
  G = G[index,index]
  y = as.vector(y[index])
  
  prefix <- "testSFI"
  
  # Run the analysis into 4 subsets and save them at a given prefix
  fm <- SFI(y,K=G,tst=1:400,trn=401:length(y),subset=c(1,4),saveAt=prefix)
  fm <- SFI(y,K=G,tst=1:400,trn=401:length(y),subset=c(2,4),saveAt=prefix)
  fm <- SFI(y,K=G,tst=1:400,trn=401:length(y),subset=c(3,4),saveAt=prefix)
  fm <- SFI(y,K=G,tst=1:400,trn=401:length(y),subset=c(4,4),saveAt=prefix)

  # Collect all results after completion
  fm <- collect(prefix)
  
  unlink(paste0(prefix,"*.RData"))   # Remove files
  unlink(paste0(prefix,"*.bin"))     # Remove files
}

\author{
Marco Lopez-Cruz (\email{lopezcru@msu.edu}) and Gustavo de los Campos
}
