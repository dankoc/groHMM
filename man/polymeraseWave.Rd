% Generated by roxygen2 (4.0.1): do not edit by hand
\name{polymeraseWave}
\alias{polymeraseWave}
\title{Given GRO-seq data, identifies the location of the polymerase 'wave' in
 up- or down- regulated genes.}
\usage{
polymeraseWave(reads1, reads2, genes, approxDist, size = 50,
  upstreamDist = 10000, TSmooth = NA, NonMap = NULL, prefix = NULL,
  emissionDistAssumption = "gamma", finterWindowSize = 10000,
  limitPCRDups = FALSE, returnVal = "simple", debug = TRUE)
}
\arguments{
\item{reads1}{Mapped reads in time point 1.}

\item{reads2}{Mapped reads in time point 2.}

\item{genes}{A set of genes in which to search for the wave.}

\item{approxDist}{The approximate position of the wave.
Suggest using 2000 [bp/ min] * time [min], for mammalian data.}

\item{size}{The size of the moving window. Suggest using 50 for direct
ligation data, and 200 for circular ligation data.  Default: 50.}

\item{upstreamDist}{The amount of upstream sequence to include
Default: 10 kb.}

\item{TSmooth}{Optimonally, outlying windows are set a maximum value
over the inter-quantile interval, specified by TSmooth.
Reasonable value: 20; Default: NA (for no smoothing).  Users are encouraged
to use this parameter ONLY in combination with the normal distribution
assumptions.}

\item{NonMap}{Optionally, un-mappable positions are trated as missing data.
NonMap passes in the list() structure for un-mappable regions.}

\item{prefix}{Optionally, writes out png images of each gene examined for
a wave.  'Prefix' denotes the file prefix for image names written to disk.
Users are encouraged to create a new directory and write in a full path.}

\item{emissionDistAssumption}{Takes values "norm", "normExp", and "gamma".
Specifies the functional form of the 'emission' distribution for states
I and II (i.e. 5' of the gene, and inside of the wave).
In our experience, "gamma" works best for highly-variable 'spikey' data,
and "norm" works for smooth data.  As a general rule of thumb, "gamma"
is used for libraries made using the direct ligation method, and "norm"
for circular ligation data.  Default: "gamma".}

\item{finterWindowSize}{Method returns 'quality' information for
each gene to which a wave was fit.  Included in these metrics are several
that define a moving window.  The moving window size is specified by
filterWindowSize.  Default: 10 kb.}

\item{limitPCRDups}{If true, counts only 1 read at each position
with >= 1 read.  NOT recommended to set this to TRUE.  Defulat: FALSE.}

\item{returnVal}{Takes value "simple" (default) or "alldata". "simple"
returns a data.frame with Pol II wave end positions.  "alldata" returns all
of the availiable data from each gene, including the full posterior
distribution of the model after EM.}

\item{debug}{If TRUE, prints error messages.}
}
\value{
Returns either a data.frame with Pol II wave end positions,
 or a List() structure with additional data, as specified by returnVal.
}
\description{
The model is a three state hidden Markov model (HMM).  States represent:
 (1) the 5' end of genes upstream of the transcription start site,
 (2) upregulated sequence, and (3) the 3' end of the gene through the
 polyadenylation site.
}
\details{
The model computes differences in read counts between the two conditions.
 Differences are assumed fit a functional form which can be specified by
 the user (using the emissionDistAssumption argument).
 Currently supported functional forms include a normal distribution
(good for GRO-seq data prepared using the circular ligation protocol),
a gamma distribution (good for 'spikey' ligation based GRO-seq data),
and a long-tailed normal+exponential distribution was implemented, but
never deployed.

Initial parameter estimates are based on initial assumptions of
 transcription rates taken from the literature.  Subsequently all parameters
 are fit using Baum-Welch expetation maximization.

Reference: Danko CG, Hah N, Luo X, Martins AL, Core L, Lis JT, Siepel A,
 Kraus WL. Signaling Pathways Differentially Affect RNA Polymerase II
 Initiation, Pausing, and Elongation Rate in Cells. Mol Cell.
 2013 Mar 19. doi:pii: S1097-2765(13)00171-8. 10.1016/j.molcel.2013.02.015.

Arguments:
}
\examples{
genes <- GRanges("chr7", IRanges(2394474,2420377), strand="+",
 SYMBOL="CYP2W1", ID="54905")
 reads1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam",
                             package="groHMM")), "GRanges")
 reads2 <- as(readGAlignments(system.file("extdata", "S40mR1.bam",
                             package="groHMM")), "GRanges")
 approxDist <- 2000*10
 # Not run:
 # pw <- polymeraseWave(reads1, reads2, genes, approxDist)
}
\author{
Charles G. Danko
}

