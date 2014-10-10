###########################################################################
##
##   Copyright 2013, 2014 Charles Danko and Minho Chae.
##
##   This program is part of the groHMM R package
##
##   groHMM is free software: you can redistribute it and/or modify it
##   under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful, but
##   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
##   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
##   for more details.
##
##   You should have received a copy of the GNU General Public License along
##   with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##########################################################################


#' writeWiggle writes a wiggle track or BigWig file suitable for uploading 
#' to the UCSC genome browser.
#'
#' @param reads GenomicRanges object representing the position of reads 
#' mapping in the genome.
#' @param file Specifies the filename for output. 
#' @param strand Takes values of "+", "-", or "*".  Computes Writes a wiggle 
#' on the speicified strand.  "*" denotes collapsing reads on both strands.  
#' Default: "*".
#' @param fileType Takes values of "wig" or "BigWig". Default: "wig".
#' @param size Size of the moving window.
#' @param normCounts A normalization factor correcting for library size 
#' or other effects.  For example, total mappible read counts might be a 
#' reasonable value.  Default: 1 (i.e. no normalization).
#' @param reverse If set to TRUE, multiplies values by -1.  
#' Used for reversing GRO-seq data on the negative (-) strand. Default: FALSE
#' @param seqinfo Seqinfo object for reads. Default: NULL.
#' @param track.type.line If set to TRUE, prints a header identifying the 
#' file as a wiggle.  Necessary to upload a custom track to the UCSC 
#' genome browser.  Default: TRUE
#' @param ...  Extra argument passed to mclapply.
#' @author Minho Chae and Charles G. Danko
#' @examples
#' S0mR1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam", 
#' package="groHMM")), "GRanges")
#' ## Not run:
#' # writeWiggle(reads=S0mR1, file="S0mR1_Plus.wig", fileType="wig", 
#' # strand="+", reverse=FALSE)
## # library(GenomicRanges)
## # si <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
## # writeWiggle(reads=S0mR1, file="S0mR1_Plus.wig", fileType="BigWig", 
## # strand="+", reverse=FALSE, seqinfo=si)
writeWiggle <- function(reads, file, strand="*", fileType="wig", size=50, 
    normCounts=NULL, reverse=FALSE, seqinfo=NULL, track.type.line=FALSE, ...) {
    W <- windowAnalysis(reads = reads, strand = strand, windowSize = size, ...)
    if (is.null(normCounts)) {
        normCounts <- 1
    }
    if (reverse) {
        normCounts <- (-1L)*normCounts
    }

    append <- FALSE
    trackWritten <- FALSE
    writeBlock <- function(chr) {
        reads_s <- GRanges(chr, 
            IRanges(seq(1, by=size, length.out=length(W[[chr]])),
            seq(size, by=size, length.out=length(W[[chr]]))),
            strand=strand, score=normCounts*as.integer(W[[chr]]))
        if (track.type.line & !trackWritten) {
            cat("track type=wiggle_0\n", file=file)
            append <<- TRUE
        }
        export(reads_s, file, format="wig", dataFormat="fixedStep", 
            append=append)
        append <<- TRUE
    }
    lapply(names(W), writeBlock)

    if (fileType=="BigWig") {
        wigToBigWig(file, seqinfo)
    }
}
