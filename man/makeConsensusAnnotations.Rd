\name{makeConsensusAnnotations}
\alias{makeConsensusAnnotations}
\title{makeConsensusAnnotations Makes a consensus annotation}
\usage{
  makeConsensusAnnotations(annotations,
    keytype = "gene_id")
}
\arguments{
  \item{annotations}{GRanges of annotations to be
  collapsed.}

  \item{keytype}{Character.  Keytype to collapse isoforms.
  Default: "gene_id"}
}
\value{
  GenomicRanges object of reads.
}
\description{
  Makes a consensus annotation which is a non overlapping
  annotation where multiply covered isoforms are used to
  represent a genomic ranges for the gene.  This also
  reduces annotations in a way to save starting site of a
  gene for overlapping gene annotations.
}
\author{
  Minho Chae
}

