\name{associateWithInterval}
\alias{associateWithInterval}
\title{associateWithInterval associates reads with the first feature that they fall inside of.}
\usage{
  associateWithInterval(features, reads)
}
\arguments{
  \item{features}{A GRanges object representing a set of
  genomic coordinates.  The meta-plot will be centered on
  the start position.}

  \item{reads}{A GRanges object representing a set of
  mapped reads.}
}
\value{
  Returns index of the feature in which a each read is
  found.  Will be a vector of integers, the same size as
  the number of reads.  NA indicates that the reads does
  not fall inside of any feature.
}
\description{
  associateWithInterval associates reads with the first
  feature that they fall inside of.
}
\author{
  Charles G. Danko and Minho Chae
}

