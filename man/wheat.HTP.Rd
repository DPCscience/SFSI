\name{wheat.HTP}
\alias{WL}
\alias{VI}
\title{High-throughput phenotypic data}
\description{

A list object type containing the hyper-spectral image data that was collected using an infrared and an hyper-spectral camera and consisted of reflectance of electromagnetic power at 250 wavelengths (or bands) within the visible and near-infrared spectrums (392-850 nm). Separate images were collected at 9 time-points covering vegetative (VEG), grain filling (GF), and maturity (MAT) stages of the crop.

The element \code{j} of the list, \code{WL[[j]]}, is the reflectance data collected at the time-point \code{j} and is a matrix with 250 columns corresponding to the wavelengths and rows matched with those of matrix \code{Y} of phenotypic records.

Likewise, the object \code{VI} is a list with 9 elements being two-columns matrices containing green and red NDVI.

Data for reflectance and GNDVI and RNDVI are reported as adjusted phenotypes from which trial, replicate and sub-block effects were removed. The resulting phenotype is the sum of a genotype and an environmental effects estimated from linear mixed models fitted within time-point.

The names (retrieved as \code{names(WL)} or \code{names(VI)}) of the list are the date (time-points) at which data was collected.

}

\source{
  International Maize and Wheat Improvement Center (CIMMYT), Mexico.
}

\keyword{datasets}