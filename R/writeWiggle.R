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


#' writeWiggle writes a wiggle track or BigWig file suitable for uploading to the UCSC genome browser.
#'
#' @param reads GenomicRanges object representing the position of reads mapping in the genome.
#' @param file Specifies the filename for output. 
#' @param strand Takes values of "+", "-", or "*".  Computes Writes a wiggle on the speicified strand.  "*" denotes collapsing reads on both strands.  Default: "*".
#' @param fileType Takes values of "wig" or "BigWig". Default: "wig".
#' @param normCounts A normalization factor correcting for library size or other effects.  For example, total mappible read counts might be a reasonable value.  Default: 1 (i.e. no normalization).
#' @param reverse If set to TRUE, multiplies values by -1.  Used for reversing GRO-seq data on the negative (-) strand. Default: FALSE
#' @param track.type.line If set to TRUE, prints a header identifying the file as a wiggle.  Necessary to upload a custom track to the UCSC genome browser.  Default: TRUE
#' @param ... Extra argument passed to export function in rtracklayer package. 
#' @return Writes a wiggle file to the specified file.
#' @author Charles G. Danko and Minho Chae
#' @examples
#' S0mR1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam", package="groHMM")), "GRanges")
#' writeWiggleNew(reads=S0mR1, file="S0mR1_Plus.wig", fileType="wig", strand="+", reverse=FALSE)
#' writeWiggleNew(reads=S0mR1, file="S0mR1_Plus.bw", fileType="BigWig", strand="+", reverse=FALSE)
writeWiggle <- function(reads, file, strand="*", fileType="wig",  
							normCounts=NULL, reverse=FALSE, track.type.line=FALSE, ...) {
	if (strand == "*") {
		reads_str <- reads
		strand(reads_str) <- "*"
		reads_str <- unique(reads_str)
	} else 
		reads_str <- unique(reads[as.character(strand(reads))==strand,])
	

	if (is.null(score(reads_str))) {
		reads_str$score <- countOverlaps(reads_str, reads) 
	}
	 
	if (reverse) 
		reads_str$score <- (-1)*reads_str$score

	if (!is.null(normCounts)) 
		reads_str$score <- normCounts*reads_str$score
	
	if (fileType=="wig") {
		if (track.type.line) {
			wigfile <- export(reads_str, format="wig", ...)
			cat("type wiggle_0\n", file=file)
			cat(wigfile, file=file, append=TRUE)
		} else 
			export(reads_str, file, format="wig", ...)
	} else 
		export(reads_str, file, format="BigWig", ...)
}
