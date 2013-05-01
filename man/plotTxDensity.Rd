\name{plotTxDensity}
\alias{plotTxDensity}
\title{plotTranscripts Plots a transcript histogram.}
\usage{
  plotTxDensity(tx, features, chrom = NA, strand = NA,
    scale = TRUE, runGenes = "best",
    brokenAnnotation = "best", filename = NA, pBar = TRUE)
}
\arguments{
  \item{tx}{GRanges of transcripts.}

  \item{features}{GRanges of non-overlapping annotatoins.}

  \item{chrom}{Character.  Target transcript chromosome; NA
  for all chromosomes. Default: NA}

  \item{strand}{Character.  "+" or "-"; NA for all
  chromosomes.  Default: NA}

  \item{scale}{Logical.  If TRUE, plot is scaled to
  relative to gene size of 30000 bp. Default: TRUE}

  \item{runGenes}{Character.  One of "best", "all", or
  "none" for the transcripts running over multiple gene
  annotations. Default: "best"}

  \item{brokenAnnotation}{Character.  One of "best", "all",
  or "none" for the transcripts breaking one annotation.
  Default: "best"}

  \item{filename}{Character.  File name for output.
  Default: NA}

  \item{pBar}{Logical.  If TRUE, progress bar is shown.
  Default: TRUE}
}
\description{
  Plots transcripts stacked together for their associated
  genes.  Gene size can be absolute or scaled to 30000 bp.
}
\author{
  Minho Chae
}

