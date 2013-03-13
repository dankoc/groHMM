\name{pausingIndex}
\alias{pausingIndex}
\title{Returns the pausing index for different genes.  TODO: DESCRIBE THE PAUSING INDEX.}
\usage{
  pausingIndex(f, p, size = 50, up = 1000, down = 1000,
    UnMAQ = NULL, debug = FALSE)
}
\arguments{
  \item{f}{data.frame of: CHR, START, END, STRAND.}

  \item{p}{data.frame of: CHR, START, END, STRAND.}

  \item{size}{The size of the moving window.}

  \item{up}{Distance upstream of each f to align and
  histogram.}

  \item{down}{Distance downstream of each f to align and
  histogram (NULL).}

  \item{UnMAQ}{Data structure representing the coordinates
  of all un-mappable regions in the genome.}

  \item{debug}{If set to TRUE, provides additional print
  options. Default: FALSE}
}
\value{
  Data.frame of the pausing indices for the input genes.
}
\description{
  Returns the pausing index for different genes.  TODO:
  DESCRIBE THE PAUSING INDEX.
}
\author{
  Charles G. Danko and Minho Chae.
}

