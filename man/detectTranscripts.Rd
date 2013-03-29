\name{detectTranscripts}
\alias{detectTranscripts}
\title{detectTranscripts detects transcripts de novo using a two-state hidden Markov model (HMM).}
\usage{
  detectTranscripts(reads = NULL, Fp = NULL, Fm = NULL,
    LtProbA = -5, LtProbB = -5, UTS = 5, size = 50,
    threshold = 0.1, debug = TRUE, ...)
}
\arguments{
  \item{reads}{A GRanges object representing a set of
  mapped reads.}

  \item{Fp}{Wiggle-formatted read counts on "+" strand.
  Optionally, Fp and Fm represent list() filled with a
  vector of counts for each chromosome.  Can detect
  transcripts starting from a fixed-step wiggle.}

  \item{Fm}{Wiggle-formatted read counts on "-" strand.}

  \item{LtProbA}{Log probability of t... .  Default: -5.
  One of these is just an initialization, and the final
  value is set by EM.  The other is a holdout parameter.}

  \item{LtProbB}{Log probability of t... .  Default: -5.}

  \item{UTS}{Varience in read counts of the untranscribed
  sequence.  Default: 5.}

  \item{size}{Log probability of t... .  Default: -5.}

  \item{threshold}{Threshold change in total likelihood,
  below which EM exits.}

  \item{debug}{If set to TRUE, provides additional print
  options. Default: FALSE}

  \item{...}{Extra argument passed to mclapply}
}
\value{
  Returns a GRanges object representing the predicted
  genomic coordinates of transcripts on both the + and -
  strand.
}
\description{
  Read counts can be specified as either a GRanges object
  (reads), or using a fixed-step wiggle-format passed in a
  list (Fp and Fm).  Either reads or BOTH Fp and Fm must be
  specified.
}
\author{
  Charles G. Danko and Minho Chae
}

