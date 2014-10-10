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

#' windowAnalysis Returns a vector of integers representing the counts of 
#' reads in a moving window.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors, set the option 'mc.cores'.
#'
#' @param reads GenomicRanges object representing the position of reads 
#' mapping in the genome.
#' @param strand Takes values of "+", "-", or "*".  "*" denotes collapsing 
#' reads on both strands.  Default: "*".
#' @param windowSize Size of the moving window. Either windowSize or 
#' stepSize must be specified.
#' @param stepSize The number of bp moved with each step.
#' @param chrom Chromosome for which to return data.  
#' Default: returns all avaliable data.
#' @param limitPCRDups Counts only one read mapping to each start site.  
#' NOTE: If set to TRUE, assumes that all reads are the same length 
#' (don't use for paired-end data).  Default: FALSE.  
#' @param ... Extra argument passed to mclapply
#' @return Returns a list object, each element of which represents a 
#' chromosome.
#' @author Charles G. Danko and Minho Chae
#' @examples
#' S0mR1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam",
#'      package="groHMM")), "GRanges")
#' ## Not run:
#' # Fp <- windowAnalysis(S0mR1, strand="+", windowSize=50)
windowAnalysis <- function(reads, strand="*", windowSize=stepSize, 
    stepSize=windowSize, chrom=NULL, limitPCRDups=FALSE, ...) {
    if (!(windowSize > 0 & (windowSize <= max(end(reads)))))
        stop("'windowSize' is out of range!")

    if (!(stepSize > 0 & (stepSize <= max(end(reads)))))
        stop("'stepSize' is out of range!")

    if (!is.null(chrom))  
        reads <- reads[seqnames(reads) == chrom,]

    readsList <- split(reads, seqnames(reads))
    if (limitPCRDups) {
        warning("Using limitPCRDups assumes all reads are the same size!  
            Don't use for paired end data!")
        readsList <- endoapply(readsList, function(x) {
            pStarts <- unique(start(x[strand(x) == "+",]))
            mEnds <- unique(end(x[strand(x) == "-",]))
            c(GRanges(seqnames(x)[1], IRanges(start=pStarts, width=1), 
                strand="+"), GRanges(seqnames(x)[1], 
                        IRanges(start=mEnds, width=1), strand="-"))

        })
    }

    # Change reads' strand 
    readsList <- endoapply(readsList, function(x) {
            if (strand == "*") 
                strand(x) <- "*"
            else
                x <- x[strand(x) == strand,]
            x 
    })

    H <- mclapply(readsList, function(x) {
            seqlevels(x) <- seqlevelsInUse(x)
            cov <- coverage(x)[[1]]
            to <- (length(cov) %/% windowSize)*windowSize
            starts <- seq(1, to, stepSize)
            vi <- Views(cov, start=starts, width=windowSize)
            Rle(viewSums(vi))
        }, ...) 
    
    
    return(H)

}
