\name{runMetaGene}
\alias{runMetaGene}
\title{Runs metagene analysis for sense and antisense direction.}
\usage{
runMetaGene(features, reads, anchorType = "TSS", size = 100L, up = 1000L,
  down = NULL, sampling = FALSE, nSampling = 1000L, samplingRatio = 0.1,
  ...)
}
\arguments{
  \item{features}{GRanges A GRanges object representing a
  set of genomic coordinates, i.e., set of genes.}

  \item{reads}{GRanges of reads.}

  \item{anchorType}{Either 'TSS' or 'TTS'.  Metagene will
  be centered on the transcription start site(TSS) or
  transcription termination site(TTS).  Default: TSS.}

  \item{size}{Numeric.  The size of the moving window.
  Default: 100L}

  \item{up}{Numeric. Distance upstream of each feature to
  align and histogram. Default: 1 kb}

  \item{down}{Numeric. Distance downstream of each feature
  to align and histogram.  If NULL, down is same as up.
  Default: NULL}

  \item{sampling}{Logical.  If TRUE, subsampling of
  Metagene is used.  Default: FALSE}

  \item{nSampling}{Numeric. Number of subsampling.
  Default: 1000L}

  \item{samplingRatio}{Numeric. Ratio of sampling for
  features.  Default: 0.1}

  \item{...}{Extra argument passed to mclapply.}
}
\value{
A integer-Rle representing the 'typical' signal centered on
a point of interest.
}
\description{
Supports parallel processing using mclapply in the
'parallel' package.  To change the number of processors set
the option 'mc.cores'.
}
\author{
Minho Chae features <- GRanges("chr7",
IRanges(start=1000:1001, width=rep(1,2)), strand=c("+",
"-")) reads <- GRanges("chr7", IRanges(start=c(1000:1003,
1100:1101), width=rep(1, 6)), strand=rep(c("+","-"), 3)) mg
<- runMetaGene(features, reads, size=4, up=10)
}

