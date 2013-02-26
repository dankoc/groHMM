#' readBed Returns GenomicRanges object after reading a file in bed format. 
#'
#' The bed files contatins seqnames, start, end, strand columns with header.
#'
#' @param file File name to read. 
#' @return GenomicRanges object of reads. 
#' @author Charles G. Danko and Minho Chae

readBed <- function(file) {
    df <- read.table(file, header=TRUE)
    return( GRanges(seqnames = Rle(df$seqnames), ranges = IRanges(df$start, df$end),
            strand = Rle(strand(df$strand))))
}

