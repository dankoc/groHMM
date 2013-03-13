\name{metaGene_nL}
\alias{metaGene_nL}
\title{Returns a histogram of the number of reads in each section of a moving window of variable size across genes.}
\usage{
  metaGene_nL(features, reads, n_windows = 1000,
    debug = FALSE)
}
\arguments{
  \item{features}{A GRanges object representing a set of
  genomic coordinates.}

  \item{reads}{A GRanges object representing a set of
  mapped reads.}

  \item{n_windows}{The number of windows to break genes
  into.}

  \item{debug}{If set to TRUE, provides additional print
  options. Default: FALSE}
}
\value{
  A vector representing the 'typical' signal across genes
  of different length.
}
\description{
  Returns a histogram of the number of reads in each
  section of a moving window of variable size across genes.
}
\author{
  Charles G. Danko and Minho Chae
}

