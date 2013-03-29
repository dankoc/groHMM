\name{metaGene}
\alias{metaGene}
\title{Returns a histogram of the number of reads in each section of a moving window centered on a certain feature.}
\usage{
  metaGene(features, reads, size, up = 1000, down = up,
    debug = FALSE, ...)
}
\arguments{
  \item{features}{A GRanges object representing a set of
  genomic coordinates.  The meta-plot will be centered on
  the start position.}

  \item{reads}{A GRanges object representing a set of
  mapped reads.}

  \item{size}{The size of the moving window.}

  \item{up}{Distance upstream of each features to align and
  histogram Default: 1 kb.}

  \item{down}{Distance downstream of each features to align
  and histogram Default: same as up.}

  \item{debug}{If set to TRUE, provides additional print
  options. Default: FALSE}

  \item{...}{Extra argument passed to mclapply}
}
\value{
  A vector representing the 'typical' signal centered on a
  point of interest.
}
\description{
  Supports parallel processing using mclapply in the
  'parallel' package.  To change the number of processors
  use the argument 'mc.cores'.
}
\author{
  Charles G. Danko and Minho Chae
}

