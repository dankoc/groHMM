\name{polymeraseWave_gamma}
\alias{polymeraseWave_gamma}
\title{Given GRO-seq data, identifies the location of the polymerase wave in up- or down-
	regulated genes.  This version is based on a full Baum-Welch EM implementation.}
\usage{
  polymeraseWave_gamma(reads1, reads2, genes, size = 50,
    approxDist, upstreamDist = 10000, TSmooth = NA,
    prefix = NULL, MinKLDiv = 1, finterWindowSize = 10000,
    debug = TRUE, returnVal = "simple")
}
\arguments{
  \item{reads1}{Mapped reads in time point 1.}

  \item{reads2}{Mapped reads in time point 2.}

  \item{genes}{A set of genes in which to search for the
  wave.}

  \item{size}{The size of the moving window. Default: 50.}

  \item{approxDist}{}

  \item{upstreamDist}{}

  \item{TSmooth}{}

  \item{prefix}{}

  \item{MinKLDiv}{}

  \item{filterWindowSize}{}

  \item{debug}{If TRUE, prints error messages.}

  \item{returnVal}{Takes value "simple" (default) or
  "alldata". "simple" returns a data.frame with Pol II wave
  end positions.  "alldata" returns all of the availiable
  data from each gene, including the full posterior
  distribution of the model after EM.}
}
\value{
  Either a data.frame with Pol II wave end positions, or a
  List() structure with additional data, as specified by
  returnVal.
}
\description{
  This is a three state HMM -- initial state representing
  the intergenic region 5' of a gene, the second
  representing the initially upregulated region, and the
  third representing the remaining sequence of a gene.
}
\details{
  We assume that upstream region is intergenic, and thus
  its emmission distriubtion is assumed to be a constant,
  set based on a representative intergenic region.  This is
  accomidated in my [1,*) framework by keeping the vairence
  constant, and scaling the mean for each gene.

  polymeraseWave_gamma assumes gamma distributed emissions
  for states two and three, and a normal distribution for
  state 1.

  Arguments:
}
\author{
  Charles G. Danko and Minho Chae
}

