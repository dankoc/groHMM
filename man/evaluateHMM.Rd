\name{evaluateHMM}
\alias{evaluateHMM}
\title{evaluateHMM Evaluates HMM calling.}
\usage{
  evaluateHMM(tx, annotations)
}
\arguments{
  \item{tx}{GRanges of transcripts predicted by HMM.}

  \item{annotations}{GRanges of non-overlapping
  annotatoins.}
}
\value{
  List of evaluation information; runGeneError,
  brokenError, overlapQuality.
}
\description{
  Evaluates HMM calling of transripts compared to know
  annotations.
}
\author{
  Minho Chae
}

