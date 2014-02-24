\name{writeWiggle}
\alias{writeWiggle}
\title{writeWiggle writes a wiggle track or BigWig file suitable for uploading to the UCSC genome browser.}
\usage{
writeWiggle(reads, file, strand = "*", fileType = "wig",
  normCounts = NULL, reverse = FALSE, track.type.line = FALSE, ...)
}
\arguments{
  \item{reads}{GenomicRanges object representing the
  position of reads mapping in the genome.}

  \item{file}{Specifies the filename for output.}

  \item{strand}{Takes values of "+", "-", or "*".  Computes
  Writes a wiggle on the speicified strand.  "*" denotes
  collapsing reads on both strands.  Default: "*".}

  \item{fileType}{Takes values of "wig" or "BigWig".
  Default: "wig".}

  \item{normCounts}{A normalization factor correcting for
  library size or other effects.  For example, total
  mappible read counts might be a reasonable value.
  Default: 1 (i.e. no normalization).}

  \item{reverse}{If set to TRUE, multiplies values by -1.
  Used for reversing GRO-seq data on the negative (-)
  strand. Default: FALSE}

  \item{track.type.line}{If set to TRUE, prints a header
  identifying the file as a wiggle.  Necessary to upload a
  custom track to the UCSC genome browser.  Default: TRUE}

  \item{...}{Extra argument passed to export function in
  rtracklayer package.}
}
\value{
Writes a wiggle file to the specified file.
}
\description{
writeWiggle writes a wiggle track or BigWig file suitable
for uploading to the UCSC genome browser.
}
\examples{
S0mR1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam", package="groHMM")), "GRanges")
writeWiggleNew(reads=S0mR1, file="S0mR1_Plus.wig", fileType="wig", strand="+", reverse=FALSE)
writeWiggleNew(reads=S0mR1, file="S0mR1_Plus.bw", fileType="BigWig", strand="+", reverse=FALSE)
}
\author{
Charles G. Danko and Minho Chae
}

