\name{associateWithInterval}
\alias{associateWithInterval}
\title{associateWithInterval associates reads with the first feature that they fall inside of.}
\usage{
associateWithInterval(features, reads, ...)
}
\arguments{
  \item{features}{A GRanges object representing a set of
  genomic coordinates.  The meta-plot will be centered on
  the start position.}

  \item{reads}{A GRanges object representing a set of
  mapped reads.}

  \item{...}{Extra argument passed to mclapply}
}
\value{
Returns index of the feature in which a each read is found.
Will be a vector of integers, the same size as the number
of reads.  NA indicates that the reads does not fall inside
of any feature.
}
\description{
Supports parallel processing using mclapply in the
'parallel' package.  To change the number of processors use
the argument 'mc.cores'.
}
\author{
Charles G. Danko and Minho Chae
}

