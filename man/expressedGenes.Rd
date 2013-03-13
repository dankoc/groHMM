\name{expressedGenes}
\alias{expressedGenes}
\title{Function identifies expressed features using the methods introduced in Core, Waterfall, Lis; Science, Dec. 2008.}
\usage{
  expressedGenes(features, reads, genomeSize = 3e+09,
    Lambda = NULL, UnMap = NULL, debug = FALSE)
}
\arguments{
  \item{features}{A GRanges object representing a set of
  genomic coordinates.  The meta-plot will be centered on
  the start position.  There can be optional "ID" column
  for gene ids.}

  \item{reads}{A GRanges object representing a set of
  mapped reads.}

  \item{genomeSize}{The size of the target genome.
  Default: 3e9, or roughly the size of the human genome.}

  \item{Lambda}{Measurement of assay noise.  Default: #
  reads/ genome size (tends to be too high for GRO-seq
  data).}

  \item{UnMap}{List object representing the position of
  un-mappable reads.  Default: not used.}

  \item{debug}{If set to true, returns the number of
  positions.  Default: FALSE.}
}
\value{
  A data.frame representing the expression p.values for
  features of interest.
}
\description{
  Function identifies expressed features using the methods
  introduced in Core, Waterfall, Lis; Science, Dec. 2008.
}
\author{
  Charles G. Danko and Minho Chae
}

