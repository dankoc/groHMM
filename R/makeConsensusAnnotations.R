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




#' makeConsensusAnnotations Makes a consensus annotation 
#'
#' Makes a non-overlapping consensus annotation.  Gene annotations are often 
#' overalpping due to #' multiple isoforms for a gene.  
#' In consensus annotation, isoforms are first reduced so that only
#' redundant intervals are used to represent a genomic interval for a gene, 
#' i.e., a gene id.
#' Remaining unresolved annotations are further reduced by truncating 3' 
#' end of annotations. 
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors, use the argument 'mc.cores'.
#'
#' @param ar GRanges of annotations to be collapsed. 
#' @param minGap Minimun gap between overlapped annotations after truncated. 
#' Default: 1L
#' @param minWidth Minimun width of consensus annotations. Default: 1000L
#' @param ... Extra argument passed to mclapply.
#' @return Returns GRanges object of annotations. 
#' @author Minho Chae
#' @examples
#' ## Not run:
#' # library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' # txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' # tx <- transcripts(txdb, vals=list(tx_chrom="chr7"), 
#' # columns=c("gene_id", "tx_id", "tx_name"))
#' # tx <- tx[grep("random", as.character(seqnames(tx)), invert=TRUE),]
#' # ca <- makeConsensusAnnotations(tx)
makeConsensusAnnotations <- function(ar, minGap=1L, minWidth=1000L, ...) {
    # check missing gene_id
    missing <- elementLengths(mcols(ar)[,"gene_id"]) == 0 
    if (any(missing)) {
        ar <- ar[!missing,]
        warning(sum(missing), " ranges do not have gene_id and they are 
            dropped")
    }

    many <- elementLengths(mcols(ar)[,"gene_id"]) > 1 
    if (any(many)) {
        ar <- ar[!many,]
        warning(sum(many), " ranges have multiple gene_id and they are 
            dropped")
    }

    ar_list <- split(ar, unlist(mcols(ar)[,"gene_id"]))
    singles <- unlist(ar_list[elementLengths(ar_list) == 1])
    isoforms <- ar_list[elementLengths(ar_list) > 1]
    
    message("Reduce isoforms(", length(isoforms),") ... ", appendLF=FALSE)
    isoforms <- GRangesList(mclapply(isoforms, function(x) {
        # For mixed strands or chrom, choose the longest    
        if ((length(seqlevelsInUse(x)) > 1) || 
                (length(unique(strand(x))) > 1)) {
            result <- x[which.max(width(x)), "gene_id"]
        } else {
            dx <- disjoin(x)
            mcols(dx)$gene_id <- mcols(x)$gene_id[1]
            olcnt <- countOverlaps(dx, x)

            multi <- dx[olcnt > 1]    # Use the disjoint ranges 
                                      # covered more than once
            if (length(multi) == 0) { # For non-overlapping isoforms, 
                                      # choose the longest
                result <- x[which.max(width(x)), "gene_id"]
            } else if (length(multi) == 1) {
                result <- multi
            } else {
                reduced <- reduce(multi)
                if (length(reduced) == 1) 
                    result <- reduced
                else (length(reduced) > 1) 
                    result <- reduced[which.max(width(reduced)),]
                
            }
            mcols(result)$gene_id <- mcols(x)$gene_id[1]
        }
        return(result)
    }, mc.cores=10))
    isoforms <- unlist(isoforms)
    message("OK")

    # Check redundancy 
    isoforms <- removeRedundant(isoforms)
    singles <- removeRedundant(singles)

    o <- findOverlaps(singles, isoforms, type="equal")
    if(length(o) != 0)
        singles <- singles[-queryHits(o),]

    o <- findOverlaps(singles, isoforms, type="within")
    if(length(o) != 0)
        singles <- singles[-queryHits(o),]

    o <- findOverlaps(isoforms, singles, type="within")
    if(length(o) != 0)
        isoforms <- isoforms[-queryHits(o),]

    noiso <- sort(c(isoforms, singles[,"gene_id"])) 
    message("Truncate overlapped ranges ... ", appendLF=FALSE)      
    # with different gene_ids
    while(!isDisjoint(noiso)) { 
        ol <- findOverlaps(noiso, ignoreSelf=TRUE, ignoreRedundant=TRUE)
        ol_gr <- GRangesList(lapply(1:length(ol), function(x) {
                        sort(c(noiso[queryHits(ol)[x]], 
                                noiso[subjectHits(ol)[x]])) 
                }))

        # Truncate 3' end
        ol_gr <- unlist(endoapply(ol_gr, function(x) {
            if (as.character(strand(x[1,])) == "+") {
                end(x[1,]) <- start(x[2,]) - minGap     
                # first range's end is truncated
            } else {
                start(x[2,]) <- end(x[1,]) + minGap     
                # sencond range's end is truncated
            }
            x
        }))

        # Remove any ranges with duplicated names since they already adujsted 
        # in the previous call
        ol_gr <- ol_gr[!duplicated(names(ol_gr)),]
        
        noiso <- noiso[-unique(c(queryHits(ol), subjectHits(ol))),] 
        # update noiso
        noiso <- c(noiso, ol_gr)
    }
    message("OK")

    noiso <- noiso[width(noiso) >= minWidth,]
    return(sort(noiso))
}

removeRedundant <- function(annox) {
    o <- findOverlaps(annox, ignoreSelf=TRUE, type="equal", 
            ignoreRedundant=TRUE)
    if(length(o) != 0)
        annox <- annox[-subjectHits(o),]

    o <- findOverlaps(annox, ignoreSelf=TRUE, type="within", 
            ignoreRedundant=TRUE)
    if(length(o) != 0)
        annox <- annox[-queryHits(o),]

    return(annox)
}

