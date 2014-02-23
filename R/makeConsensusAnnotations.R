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





##########################################################################
##      makeConsensusAnnotations
##      Date: 2014-2-19 
#' makeConsensusAnnotations Makes a consensus annotation 
#'
#' Makes a non-overlapping consensus annotation.  Gene annotations are often overalpping due to 
#' multiple isoforms for a gene.  In consensus annotation, isoforms are first reduced so that only
#' redundant intervals are used to represent a genomic interval for a gene, i.e., a gene id.
#' Remaining unresolved annotations are further reduced by truncating 3' end of annotations. 
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  To change the number of processors
#' use the argument 'mc.cores'.
#'
#' @param ar GRanges of annotations to be collapsed. 
#' @param minGap Minimun gap between overlapped annotations after truncated. Default: 1L
#' @return GenomicRanges object of annotations. 
#' @author Minho Chae
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' tx <- transcripts(txdb, vals=list(tx_chrom="chr7"), columns=c("gene_id", "tx_id", "tx_name"))
#' tx <- tx[grep("random", as.character(seqnames(tx)), invert=TRUE),]
#' ca <- makeConsensusAnnotations(tx)
##########################################################################
makeConsensusAnnotations <- function(ar, minGap=1L, ...) {
	# check missing gene_id
	missing <- elementLengths(mcols(ar)[,"gene_id"]) == 0 
	if (any(missing)) {
		ar <- ar[!missing,]
		warning(sum(missing), " ranges do not have gene_id and they are dropped")
	}

	many <-	elementLengths(mcols(ar)[,"gene_id"]) > 1 
	if (any(many)) {
		ar <- ar[!many,]
		warning(sum(many), " ranges have multiple gene_id and they are dropped")
	}

	ar_list <- split(ar, unlist(mcols(ar)[,"gene_id"]))
	singles <- unlist(ar_list[elementLengths(ar_list) == 1])
	isoforms <- ar_list[elementLengths(ar_list) > 1]
	
	message("Reduce isoforms(", length(isoforms),") ... ", appendLF=FALSE)

	noiso <- GRangesList(mclapply(isoforms, function(x) { 
		dx <- disjoin(x)
		mcols(dx)$gene_id <- mcols(x)$gene_id[1]
		olcnt <- countOverlaps(dx, x) 
			
		multi_dx <- dx[olcnt > 1] 		# Use disjoint ranges covered more than once
		if (length(multi_dx) == 0) {
			multi_dx <- dx
		} else if (length(multi_dx) == 1) {
			return(multi_dx)
		} else {	
			# In order to preserve TSS of the annotations, 
			# choose 5' or 3' most one for plus or minus strand, respectively.
			multi_dx_str <- unique(as.character(strand(multi_dx))) 
			if (length(multi_dx_str) == 1) {  	
				multi_dx <- sort(multi_dx)      
											
				if (multi_dx_str == "+") {
					return(multi_dx[1,])
				} else {
					return(multi_dx[length(multi_dx),])
				}	
			} else {  # for mixed strands, choose the longest
				return(multi_dx[which.max(width(multi_dx)),])
			}
		}
	},...))
	noiso <- unlist(noiso)
	message("OK")

	noiso <- sort(c(noiso, singles[,"gene_id"])) 
	# query within or equal to subject "
	ol_w <- findOverlaps(noiso, type="within", ignoreSelf=TRUE, ignoreRedundant=TRUE) 
	if (length(ol_w)) {
		noiso <- noiso[-unique(queryHits(ol_w)),]
	}
	
	message("Truncate overlapped ranges ... ", appendLF=FALSE)		# with different gene_ids
	while(!isDisjoint(noiso)) { 
		ol <- findOverlaps(noiso, ignoreSelf=TRUE, ignoreRedundant=TRUE)
		ol_gr <- GRangesList(lapply(1:length(ol), function(x) {
						sort(c(noiso[queryHits(ol)[x]], noiso[subjectHits(ol)[x]])) 
					}))

		ol_gr <- unlist(endoapply(ol_gr, function(x) {
			if (as.character(strand(x[1])) == "+") {
				end(x[1]) <- start(x[2]) - minGap 	# first range's end is truncated
			} else {
				start(x[2]) <- end(x[1]) + minGap 	# sencond range's end is truncated
			}
			x
		}))

		# Remove any ranges with duplicated names since they already adujsted in the previous call
		ol_gr <- ol_gr[!duplicated(names(ol_gr)),]
		
		noiso <- noiso[-unique(c(queryHits(ol), subjectHits(ol))),] # update noiso
		noiso <- c(noiso, ol_gr)
	}
	message("OK")

	return(sort(noiso))


}
	#H <- mclapply(chrom, makeConsensusAnnotations_foreachChrom, ar=ar, minGap=minGap, ...)
	#return(sort(unlist(GRangesList(H))))
