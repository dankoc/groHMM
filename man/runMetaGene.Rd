\name{runMetaGene}
\alias{runMetaGene}
\title{runMetaGene Runs Metagene analysis.}
\usage{
  runMetaGene(features, reads, size = 100, up = 10000,
    down = NULL, normCounts = 1, sampling = FALSE,
    nSampling = 1000, debug = FALSE)
}
\arguments{
  \item{features}{GRanges whose length is 1, i.e.,
  Transcription Start Site (TSS)}

  \item{reads}{GRanges of reads.}

  \item{size}{Numeric.  The size of the moving window.
  Default: 100}

  \item{up}{Numeric. Distance upstream of each fgr to align
  and histogram. Default: 10000}

  \item{down}{Numeric. Distance downstream of each fgr to
  align and histogram.  If NULL, down is same as up.
  Default: NULL}

  \item{normCounts}{Numeric.  Reads are multiplied by
  normCounts.Thhe size of the moving window. Default: 1}

  \item{sampling}{Logical.  If TRUE, subsampling of
  Metagene is used.  Default: FALSE}

  \item{nSampling}{Numeric. Number of subsampling.
  Default: 1000}

  \item{debug}{Logical. If set to TRUE, provides additional
  print options. Default: FALSE}
}
\value{
  List of vectors representing the 'typical' signal
  centered on the genomic features of interest.
}
\description{
  Returns a histogram of the number of reads in each
  section of a moving window centered on a certain feature
  for both sense and anntisense directions.  It has option
  for subsampling.
}
\author{
  Minho Chae
}

