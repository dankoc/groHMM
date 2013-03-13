\name{limitToXkb}
\alias{limitToXkb}
\title{limitToXkb truncates a set of genomic itnervals at a constant, maximum size.}
\usage{
  limitToXkb(features, offset = 1000, size = 13000)
}
\arguments{
  \item{features}{A GRanges object representing a set of
  genomic coordinates.  The meta-plot will be centered on
  the start position.}

  \item{offset}{Starts the interval from this position
  relative to the start of each genomic features.}

  \item{size}{Specifies the size of the window.}
}
\value{
  Returns a new 'GRanges' object representing the size new
  size.
}
\description{
  limitToXkb truncates a set of genomic itnervals at a
  constant, maximum size.
}
\author{
  Minho Chae and Charles G. Danko
}

