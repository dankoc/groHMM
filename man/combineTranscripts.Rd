\name{combineTranscripts}
\alias{combineTranscripts}
\title{combineTranscripts Combines transnscipts.}
\usage{
  combineTranscripts(tx, annotations, geneSize = 1000,
    threshold = 0.8, debug = FALSE)
}
\arguments{
  \item{tx}{GRanges of transcripts.}

  \item{annotations}{GRanges of non-overlapping annotations
  for reference.}

  \item{geneSize}{Numeric. Minimum gene size in annotations
  to be used as reference. Default: 1000}

  \item{threshold}{Numeric. Threshold for calling the gene
  part of the transcript.  Default: 0.8}

  \item{debug}{Logical.  If set to TRUE, show easch step in
  a plot. Default: FALSE}
}
\description{
  Combines transcripts that are within the same gene
  annotation, combining smaller transcripts for genes with
  low regulation into a single transcript representing the
  gene.
}
\author{
  Minho Chae and Charles G. Danko
}

