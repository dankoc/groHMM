\name{plotTranscripts}
\alias{plotTranscripts}
\title{plotTranscripts Plots transcripts.}
\usage{
  plotTranscripts(tx, features, chrom = NA, strand = NA,
    scale = TRUE, runGenes = "best",
    brokenAnnotation = "best", first = NA, last = NA,
    filename = NA)
}
\arguments{
  \item{tx}{GRanges of transcripts.}

  \item{features}{GRanges of non-overlapping annotatoins.}

  \item{chrom}{Character.  Target transcript chromosome; NA
  for all chromosomes.  Default: NA}

  \item{strand}{Character.  "+" or "-"; NA for all
  chromosomes.  Default: NA}

  \item{scale}{Logical.  If TRUE, plot is scaled to
  relative to gene size of 30000 bp. Default: TRUE}

  \item{runGenes}{Character.  One of "best", "all", or
  "none" for the transcripts running over multiple gene
  annotations. Default: "best"}

  \item{brokenAnnotation}{Character.  One of "best", "all",
  or "none" for the transcripts breaking one annotation.
  Defualt: "best"}

  \item{first}{Integer. First n transcripts to plot after
  transcripts are ordered by their gene sizes. Default: NA}

  \item{last}{Integer. Last n transcripts to plot after
  transcripts are ordered by their gene sizes. Default: NA}

  \item{filename}{Character.  File name for output.
  Default: NA}
}
\description{
  Plots transcripts with their associated genes.
}
\author{
  Minho Chae
}

