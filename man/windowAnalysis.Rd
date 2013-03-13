\name{windowAnalysis}
\alias{windowAnalysis}
\title{windowAnalysis Returns a vector of integers representing the counts of reads in a moving window.}
\usage{
  windowAnalysis(reads, strand = "N",
    window_size = (step_size - 1),
    step_size = (window_size + 1), chrom = NULL, start = 0,
    end = NULL, limitPCRDups = FALSE, debug = FALSE)
}
\arguments{
  \item{reads}{GenomicRanges object representing the
  position of reads mapping in the genome.}

  \item{strand}{Takes values of "+", "-", or "N".  Computes
  Writes a wiggle on the speicified strand.  "N" denotes
  collapsing reads on both strands.  Default: "N".}

  \item{window_size}{Size of the moving window. Either
  window_size or step_size must be specified.}

  \item{step_size}{The number of bp moved with each step.}

  \item{chrom}{Chromosome for which to return data.
  Default: returns all avaliable data.}

  \item{start}{The start position in genomic coordinates.
  Default: returns all avaliable data.}

  \item{end}{The end position in genomic coordinates.
  Default: returns all avaliable data.}

  \item{limitPCRDups}{Counts only one read mapping to each
  start site.  NOTE: If set to TRUE, assumes that all reads
  are the same length (don't use for paired-end data).
  Default: FALSE.}

  \item{debug}{If set to TRUE, provides additional print
  options. Default: FALSE}
}
\value{
  List object, each element of which represents a
  chromosome.
}
\description{
  windowAnalysis Returns a vector of integers representing
  the counts of reads in a moving window.
}
\author{
  Charles G. Danko and Minho Chae
}

