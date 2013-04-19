\name{writeWiggle}
\alias{writeWiggle}
\title{writeWiggle writes a wiggle track suitable for uploading to the UCSC genome browser.
Currently only supports writing a fixed-step wiggle.}
\usage{
  writeWiggle(reads, file, strand = "N", size = 50,
    normCounts = 1, sep.chrom = FALSE, reverse = FALSE,
    track.type.line = TRUE, debug = FALSE, ...)
}
\arguments{
  \item{reads}{GenomicRanges object representing the
  position of reads mapping in the genome.}

  \item{file}{Specifies the file prefix for the output
  wiggle.}

  \item{strand}{Takes values of "+", "-", or "N".  Computes
  Writes a wiggle on the speicified strand.  "N" denotes
  collapsing reads on both strands.  Default: "N".}

  \item{size}{Size of non-overlapping windows to write.
  Default: 50 bp.}

  \item{normCounts}{A normalization factor correcting for
  library size or other effects.  For example, total
  mappible read counts might be a reasonable value.
  Default: 1 (i.e. no normalization).}

  \item{sep.chrom}{If set to TRUE, will write a separate
  wiggle file for each chromosome.  Default: FALSE}

  \item{reverse}{If set to TRUE, multiplies values by -1.
  Used for reversing GRO-seq data on the negative (-)
  strand. Default: FALSE}

  \item{track.type.line}{If set to TRUE, prints a header
  identifying the file as a wiggle.  Necessary to upload a
  custom track to the UCSC genome browser.  Default: TRUE}

  \item{debug}{If set to TRUE, provides additional print
  options. Default: FALSE}

  \item{...}{Extra argument passed to mclapply}
}
\value{
  Writes a wiggle file to the specified file.
}
\description{
  writeWiggle writes a wiggle track suitable for uploading
  to the UCSC genome browser. Currently only supports
  writing a fixed-step wiggle.
}
\author{
  Charles G. Danko
}

