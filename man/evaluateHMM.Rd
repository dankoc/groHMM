\name{evaluateHMM}
\alias{evaluateHMM}
\title{evaluateHMM Evaluates HMM calling.}
\usage{
  evaluateHMM(tr, annotations)
}
\arguments{
  \item{tr}{GRanges of transcripts predicted by HMM.}

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

