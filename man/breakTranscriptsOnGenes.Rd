\name{breakTranscriptsOnGenes}
\alias{breakTranscriptsOnGenes}
\title{breakTranscriptsOnGenes Breaks transcripts on genes}
\usage{
  breakTranscriptsOnGenes(tx, annotations, strand = "+",
    geneSize = 5000, threshold = 0.8, gap = 5,
    debug = FALSE)
}
\arguments{
  \item{tx}{GRanges of transcripts.}

  \item{annotations}{GRanges of non-overlapping annotations
  for reference.}

  \item{strand}{Takes "+" or "-" Default: "+"}

  \item{geneSize}{Numeric. Minimum gene size in gr to be
  used as reference. Default: 5000}

  \item{threshold}{Numeric. Threshold for calling the gene
  part of the transcript.  Default: 0.8}

  \item{gap}{Numeric.  Gap (bp) between broken transcripts.
  Default: 5}

  \item{debug}{Logical.  If set to TRUE, show easch step in
  a plot. Default: FALSE}
}
\description{
  Breaks transcripts when they are fully overlap multiple
  well annotated genes.
}
\author{
  Minho Chae and Charles G. Danko
}

