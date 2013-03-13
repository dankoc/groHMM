\name{countMappableReadsInInterval}
\alias{countMappableReadsInInterval}
\title{countMappableReadsInInterval counts the number of mappable reads in a set of genomic features.}
\usage{
  countMappableReadsInInterval(features, UnMap,
    debug = FALSE)
}
\arguments{
  \item{features}{A GRanges object representing a set of
  genomic coordinates.  The meta-plot will be centered on
  the start position.}

  \item{UnMap}{List object representing the position of
  un-mappable reads.  Default: not used.}

  \item{debug}{If set to TRUE, provides additional print
  options. Default: FALSE}
}
\value{
  Returns a vector of counts, each representing the number
  of reads inside each genomic interval.
}
\description{
  countMappableReadsInInterval counts the number of
  mappable reads in a set of genomic features.
}
\author{
  Charles G. Danko and Minho Chae
}

