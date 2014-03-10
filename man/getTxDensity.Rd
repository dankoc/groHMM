\name{getTxDensity}
\alias{getTxDensity}
\title{getTxDensity Calculates transcript density.}
\usage{
getTxDensity(tx, annox, plot = TRUE, scale = 1000L, nSampling = 0L,
  samplingRatio = 0.1, ...)
}
\arguments{
  \item{tx}{GRanges of transcripts.}

  \item{annox}{GRanges of non-overlapping annotatoins.}

  \item{plot}{Logical.  If TRUE, plot transcript density.
  Default: TRUE}

  \item{strand}{Character.  "+" or "-"; NA for all
  chromosomes.  Default: NA}

  \item{scale}{Numeric. Scaled size of a gene for
  transcript density calculation. Default: 1000L}

  \item{nSampling}{Numeric. Number of subsampling.
  Default: 0L}

  \item{samplingRatio}{Numeric. Ratio of sampling for
  annotations.  Default: 0.1}

  \item{...}{Extra argument passed to mclapply.}
}
\value{
A list of FTD, TTD, PostTTS, and AUC.
}
\description{
Calculates transcript density for transcripts which
overlapps with annotations. For 'run genes together' or
'broken up a single annotation' errors, best overlapped
transcripts or annotations are used.
}
\examples{
tx <- GRanges("chr7", IRanges(start=seq(1000,4000, by=1000), width=seq(1000, 1300, by=100)), strand=rep("+", 4))
annox <- GRanges("chr7", IRanges(start=seq(1100,4100, by=1000), width=seq(900, 1200, by=100)), strand=rep("+", 4))
density <- plotTxDensity(tx, annox)
}
\author{
Minho Chae
}

