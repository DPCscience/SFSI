\name{wheat}
\alias{wheat.E1}
\alias{wheat.E2}
\alias{wheat.E3}
\alias{wheat.E4}
\alias{wheat.G}
\title{Wheat dataset}
\description{

The data set consists of 1,092 inbred wheat lines grouped into 39 trials and grown during the 2013-2014 season at the Norman Borlaug experimental research station in Ciudad Obregon, Sonora, Mexico. Each trial consisted of 28 breeding lines that were arranged in an alpha-lattice design with three replicates and six sub-blocks. The trials were grown in four different environments: 
\itemize{
\item E1: Flat-Drought (sowing in flat with irrigation of 180 mm through drip system)
\item E2: Bed-2IR (sowing in bed with 2 irrigations approximately 250 mm)
\item E3: Bed-5IR (bed sowing with 5 normal irrigations)
\item E4: Bed-EHeat (bed sowing 30 days before optimal planting date with 5 normal irrigations approximately 500 mm)
}

Reflectance phenotypic data were collected from the fields using both infrared (A600 series Infrared camera, FLIR, Wilsonville, OR) and hyper-spectral (A-series, Micro-Hyperspec, VNIR Headwall Photonics, Fitchburg, MA) cameras mounted on a Piper PA-16 Clipper aircraft on 9 different dates (time-points) between January 10 and March 27th, 2014. During each flight, data from 250 wavelengths ranging from 392 to 850 nm were collected for each pixel in the pictures. The average reflectance of all the pixels for each wavelength was calculated from each of the geo-referenced trial plots and reported as each line reflectance.

Data was used in Lopez-Cruz et al. (2019)
}

\usage{
  data(wheat.E1)
  data(wheat.E2)
  data(wheat.E3)
  data(wheat.E4)
  data(wheat.G)
}

\format{
Whitin each dataset, matrix Y contains the phenotypic data for grain yield (YLD), days to heading (DTH), days to maturity (DTM), and plant height (PH). The object WL is a list that contains for each time-point the hyperspectral reflectance data. Object VI is a list that contains for each time-point green and red NDVI obtained from reflectance data. The matrix G contains genomic relationships for a subset of 807 lines computed from SNP markers. 
}

\source{
  International Maize and Wheat Improvement Center (CIMMYT), Mexico.
}

\references{
Lopez-Cruz M, Olson E, Rovere G, Crossa J, Dreisigacker S, Suchismita M, Singh R, de los Campos G (2019). Regularized selection indices for breeding value prediction using hyper-spectral image data. \emph{preprint bioRxiv}.
}
\keyword{datasets}