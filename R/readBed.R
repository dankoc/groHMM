#' readBed Returns a GenomicRanges object constrcuted from the specified bed file. 
#'
#' Bed file format is assumed to be either four column: seqnames, start, end, strand columns; or six column: seqnames, start, end, name, score, and strand.  
#' There is currently no support for bed files without strand information.  
#'
#' Any additional arguments availiable to read.table can be specified.  
#'
#' @param file Path to the input file.
#' @return GenomicRanges object, representing mapped reads. 
#' @author Minho Chae and Charles G. Danko.
readBed <- function(file, ...) {
    df <- read.table(file, ...)
	## Do we want an option for a 3-column bed (i.e. strand not specified)?  It's technically a legal UCSC format.  
	if(NCOL(df) == 4) col.names(df) <- c("seqnames", "start", "end", "strand")
    if(NCOL(df) == 6) col.names(df) <- c("seqnames", "start", "end", "name", "score", "strand")
    return( GRanges(seqnames = Rle(df$seqnames), ranges = IRanges(df$start, df$end),
            strand = Rle(strand(df$strand))))
}

