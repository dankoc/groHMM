###########################################################################
##
##   Copyright 2012, 2013 Minho Chae.
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


#' readBed Returns a GenomicRanges object constrcuted from the specified bed file. 
#'
#' Bed file format is assumed to be either four column: seqnames, start, end, strand columns; or six column: seqnames, start, end, name, score, and strand.  
#' Three column format is also possible when there is no strand information. 
#'
#' Any additional arguments availiable to read.table can be specified.  
#'
#' @param file Path to the input file.
#' @param ... Extra argument passed to read.table 
#' @return GenomicRanges object, representing mapped reads. 
#' @author Minho Chae and Charles G. Danko.
readBed <- function(file, ...) {
    df <- read.table(file, ...)
	if(NCOL(df) == 3) {
		colnames(df) <- c("seqnames", "start", "end")
		df <- cbind(df, strand=Rle("*", NROW(df)))
	}
	if(NCOL(df) == 4) colnames(df) <- c("seqnames", "start", "end", "strand")
	if(NCOL(df) == 6) colnames(df) <- c("seqnames", "start", "end", "name", "score", "strand")
	return( GRanges(seqnames = Rle(df$seqnames), ranges = IRanges(df$start, df$end),
            strand = Rle(strand(df$strand))))
}

