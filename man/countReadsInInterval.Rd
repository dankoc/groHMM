\name{countReadsInInterval}
\alias{countReadsInInterval}
\title{countReadsInInterval counts the number of reads that fall inside of each genomic feature.}
\usage{
  countReadsInInterval(features, reads)
}
\arguments{
  \item{features}{A GRanges object representing a set of
  genomic coordinates.  The meta-plot will be centered on
  the start position.}

  \item{reads}{A GRanges object representing a set of
  mapped reads.}
}
\value{
  Returns a vector of counts, each representing the number
  of reads inside each genomic interval.
}
\description{
  Only counts reads on the appropriate strand.
}
\author{
  Charles G. Danko and Minho Chae
}

