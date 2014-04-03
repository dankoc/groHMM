\name{groHMM-package}
\alias{groHMM-package}
\alias{groHMM}
\docType{package}
\title{
groHMM: GRO-seq Analysis Pipeline
}
\description{
groHMM was developed for analysis of GRO-seq data, which provides a genome wide 'map' of the position and orientation of all
transcriptionally active RNA polymerases. groHMM predicts the boundaries of
transcriptional activity across the genome de novo using a two-state hidden Markov
model (HMM).  The model essentially divides the genome into 'transcribed' and
'non-transcribed' regions in a strand specific manner. 

We also use HMMs to identify the leading edge of Pol II at genes activated by a
stimulus in GRO-seq time course data.  This approach allows the genome-wide
interrogation of transcription rates in cells.

In addition to these advanced features, groHMM provides wrapper functions for
counting raw reads, generating wiggle files for visualization,
and creating metagene (averaging) plots.  Although groHMM is tailored
towards GRO-seq data, the same functions and analytical methodologies can, in principal,
be applied to a wide variety of other short read data sets.
}
\details{
\tabular{ll}{
Package: \tab groHMM\cr
Type: \tab Package\cr
Version: \tab 0.99.0\cr
Date: \tab 2014-04-02\cr
License: \tab GPL (>=3)\cr
LazyLoad: \tab yes\cr
Depends: \tab R (>= 2.14.0), MASS, GenomicRanges, rtracklayer, parallel\cr
}
}
\author{
Charles G. Danko, Minho Chae, Andre Martins  \cr

Maintainer: Minho Chae<minho.chae@gmail.com>
}
\references{
Luo, X., Chae, M., Krishnakumar, R,, Danko, C., Kraus, L. Dynamic reorganization of the AC16 cardiomyocyte transcriptome in response to TNFa
        signaling revealed by integrated genomic analyses. BMC Genomics. 2014 Feb 24;15(1):155

Hah, N., Danko, C., Core, L., Waterfall, J., Siepel, A., Lis, J., Kraus, L.  A Rapid, Extensive, and Transient Transcriptional 
Response to Estrogen Signaling in Breast Cancer Cells. Cell. 2011 May 13;145(4):622-34 
}
\keyword{ package }
