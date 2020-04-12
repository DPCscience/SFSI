\name{G}
\alias{G}
\title{Genomic relationship matrix}
\description{

All 1,092 lines, along with a larger set of 19,965 lines, were sequenced using GBS (genotyping-by-sequencing) technology at 192-plexing on Illumina HiSeq2000 or HiSeq2500 with 1 x 100 bp reads. Single nucleotide polymorphisms (SNP) were called across all lines anchored to the genome assembly of Chinese Spring (International Wheat Genome Sequencing Consortium 2014). 

Markers were filtered so that the percent of missing data and percent heterozygosity per marker did not exceed 80\% and 20\%, respectively. Next, SNP were extracted for the 1,092 lines and then filtered so that lines >50\% missing data were removed. Markers were recoded as –1, 0, and 1, corresponding to homozygous for the minor allele, heterozygous, and homozygous for the major allele respectively. Next, markers with a minor allele frequency >0.01 and >50\% of missing data were removed. Remaining markers with missing data were imputed with the marker mean. 

A total of 2,334 SNPs scored on 807 individuals was used to obtain the genomic relationship matrix acording to VanRaden (2008).
}

\source{
  International Maize and Wheat Improvement Center (CIMMYT), Mexico.
}

\references{
VanRaden PM (2008). Efficient methods to compute genomic predictions. \emph{Journal of Dairy Science}, \bold{91}(11), 4414–4423.
}
\keyword{datasets}