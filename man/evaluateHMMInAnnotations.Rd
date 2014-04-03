\name{evaluateHMMInAnnotations}
\alias{evaluateHMMInAnnotations}
\title{evaluateHMM Evaluates HMM calling.}
\usage{
evaluateHMMInAnnotations(tx, annox)
}
\arguments{
  \item{tx}{GRanges of transcripts predicted by HMM.}

  \item{annox}{GRanges of non-overlapping annotatoins.}
}
\value{
a list of error information; runGenes,
brokenSingleAnnotations, total, and rate.
}
\description{
Evaluates HMM calling of transripts compared to known
annotations.
}
\examples{
tx <- GRanges("chr7", IRanges(start=seq(100, 1000, by=200), width=seq(100, 1000, by=100)), strand="+")
annox <- GRanges("chr7", IRanges(start=seq(110, 1100, by=150), width=seq(100, 1000, by=150)), strand="+")
error <- evaluateHMMAnnotations(tx, annox)
}
\author{
Minho Chae
}

