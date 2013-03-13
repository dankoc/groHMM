\name{metaGeneMatrix}
\alias{metaGeneMatrix}
\title{Returns a matrix, with rows representing read counts across a specified gene, or other features of interest.}
\usage{
  metaGeneMatrix(features, reads, size = 50, up = 1000,
    down = up, debug = FALSE)
}
\arguments{
  \item{features}{A GRanges object representing a set of
  genomic coordinates.}

  \item{reads}{A GRanges object representing a set of
  mapped reads.}

  \item{size}{The size of the moving window.}

  \item{up}{Distance upstream of each f to align and
  histogram Default: 1 kb.}

  \item{down}{Distance downstream of each f to align and
  histogram Default: same as up.}

  \item{debug}{If set to TRUE, provides additional print
  options. Default: FALSE}
}
\value{
  A vector representing the 'typical' signal across genes
  of different length.
}
\description{
  Returns a matrix, with rows representing read counts
  across a specified gene, or other features of interest.
}
\author{
  Charles G. Danko and Minho Chae
}

