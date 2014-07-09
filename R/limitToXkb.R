###########################################################################
##
##   Copyright 2009, 2010, 2011, 2012, 2013 Charles Danko and Minho Chae.
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

#' limitToXkb truncates a set of genomic itnervals at a constant, maximum size.
#'
#' @param features A GRanges object representing a set of genomic coordinates. 
#' The meta-plot will be centered on the start position.
#' @param offset Starts the interval from this position relative to the start 
#' of each genomic features.
#' @param size Specifies the size of the window.
#' @return Returns GRanges object with new genomic coordiates.
#' @author Minho Chae and Charles G. Danko
#' @examples
#' tx <- GRanges("chr7", IRanges(1000, 30000), strand="+")
#' newTX <- limitToXkb(tx)
##  This function limits a genomic range to a samll region relative to the 
## transcription site.
limitToXkb <- function(features, offset=1000, size=13000) {
    w <- width(features)

    # 1. do nothing for w < offset 
    # 2. offset < w and w < size
    small  <- (offset < w) & (w < size)
    if (any(small)) { 
        features[small,] <- flank(features[small,], -1*(w[small]-offset), 
            start=FALSE)
    }

    # 2. w > size
    big  <- w > size 
    if (any(big)) {
        features[big,] <- resize(features[big,], width=size)

        bigPlus <- big & as.character(strand(features))=="+"
        if (any(bigPlus)) 
            start(features[bigPlus,]) <- start(features[bigPlus,]) + offset 

        bigMinus <- big & as.character(strand(features))=="-"
        if (any(bigMinus)) 
            end(features[bigMinus,]) <- end(features[bigMinus,]) - offset 
    }

    return(features)
}

#' countMappableReadsInInterval counts the number of mappable reads in a set 
#' of genomic features.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors, use the argument 'mc.cores'.
#'
#' @param features A GRanges object representing a set of genomic coordinates.
#' The meta-plot will be centered on the start position.
#' @param UnMap List object representing the position of un-mappable reads.  
#' Default: not used.
#' @param debug If set to TRUE, provides additional print options. 
#' Default: FALSE
#' @param ... Extra argument passed to mclapply
#' @return Returns a vector of counts, each representing the number of reads 
#' inside each genomic interval.
#' @author Charles G. Danko and Minho Chae
##
##  This function takes information from BED file to represent regions 
## (as in CountReadsInInterval), and  a list structure representing unmappable 
## positions.  Counts the number of mappable positions in  the interval.
##      f  == features in a BED format, where columns represent: 
##            Chr, Start, End, Strand, ID.
##      UnMAQ  == List structure of the un-mappable bases in the genome.
##
##  Addumptions:
##  (1) FeatureStart and FeatureEnd should be vectors ONLY on the same 
##  chromosome as PROBEStart!!!
##  (2) Take care: only returns probes on the same strand as the feature, f.
##  to return all probes, force all strands to be the same.
countMappableReadsInInterval <- function(features, UnMap, debug=FALSE, ...) {

    C <- sort(unique(as.character(seqnames(features))))
    
    ## Run parallel version.
    mcp <- mclapply(seq_along(C), countMappableReadsInInterval_foreachChrom, 
                    C=C, features=features, UnMap=UnMap, ...)

    ## Convert to a vector.
    FT <- rep(0,NROW(features))
    for(i in seq_along(C)) {
        indxF   <- which(as.character(seqnames(features)) == C[i])

        if(NROW(indxF) >0) {
            FT[indxF][mcp[[i]][["ord"]]] <- as.integer(mcp[[i]][["Difference"]])
        }
    }

    return(FT)
}


countMappableReadsInInterval_foreachChrom <- function(i, C, features, UnMap) {
    indxF   <- which(as.character(seqnames(features)) == C[i])

    if(NROW(indxF) >0) {
        # Order -- Make sure, b/c this is one of our main assumptions.  
        # Otherwise violated for DBTSS.
        Ford <- order(start(features[indxF,]))
        
        # Type coersions.
        FeatureStart    <- start(features[indxF,][Ford])
        FeatureEnd  <- end(features[indxF,][Ford])
        FeatureStr  <- as.character(strand(features[indxF,][Ford]))

        # Set dimensions.
        dim(FeatureStart)   <- c(NROW(FeatureStart), NCOL(FeatureStart))
        dim(FeatureEnd)     <- c(NROW(FeatureEnd),   NCOL(FeatureEnd))
        dim(FeatureStr)     <- c(NROW(FeatureStr),   NCOL(FeatureStr))
        
        ## Count start index.
        chr_indx <- which(UnMap[[1]][[1]] == C[i])
        CHRSIZE <- as.integer(UnMap[[1]][[2]][chr_indx])
        CHRSTART <- as.integer(0)
        if(chr_indx > 1) {  ## Running on 1:0 gives c(1, 0)
            CHRSTART <- as.integer( 
                sum(UnMap[[1]][[2]][
                    c(1:(chr_indx-1))
                ]) +1)
        }

        if(debug) {
            message(C[i],": Counting unMAQable regions.")
            message("CHRSIZE:", CHRSIZE, "CHRSTART:", CHRSTART)
        }

        ## Count unMAQable regions, and size of everything ... 
        nonmappable <- .Call("CountUnMAQableReads", FeatureStart, FeatureEnd, 
                UnMap[[2]], CHRSTART, CHRSIZE, PACKAGE = "groHMM")

        ## Adjust size of gene body.
        Difference <- (FeatureEnd - FeatureStart) - nonmappable + 1 
        ## Otherwise, get -1 for some.

        if(debug) {
            print(head(nonmappable))
            print(as.integer(head(Difference)))
        }

        #F[indxF][Ford] <- as.integer(Difference)
        return(list(Difference= Difference, ord= Ford))
    }
    return(integer(0))
}



 
#' readBed Returns a GenomicRanges object constrcuted from the specified bed 
#' file.
#'
#' Bed file format is assumed to be either four column: seqnames, start, end, 
#' strand columns; or six column: seqnames, start, end, name, score, and strand.
#' Three column format is also possible when there is no strand information.
#'
#' Any additional arguments availiable to read.table can be specified.
#'
#' @param file Path to the input file.
#' @param ... Extra argument passed to read.table
#' @return Returns GRanges object representing mapped reads.
#' @author Minho Chae and Charles G. Danko.
readBed <- function(file, ...) {
    df <- read.table(file, ...)
        if(NCOL(df) == 3) {
                colnames(df) <- c("seqnames", "start", "end")
                df <- cbind(df, strand=Rle("*", NROW(df)))
        }
        if(NCOL(df) == 4) colnames(df) <- c("seqnames", "start", "end", 
                                            "strand")
        if(NCOL(df) == 6) colnames(df) <- c("seqnames", "start", "end", 
                                            "name", "score", "strand")
        return( GRanges(seqnames = Rle(df$seqnames), ranges = 
                IRanges(df$start, df$end), strand = Rle(strand(df$strand))))
}
