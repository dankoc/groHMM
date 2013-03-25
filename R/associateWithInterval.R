###########################################################################
##
##   Copyright 2009, 2010, 2011 Charles Danko.
##
##   This program is part of the GRO-seq R package
##
##   GRO-seq is free software: you can redistribute it and/or modify it 
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


########################################################################
##
##	FeatureStart and FeatureEnd should be vectors ONLY on the same 
##	chromosome as PROBEStart!!!
##
##	--Calls C function for speed.
##	--Returns index of feature island to which a probe belongs, assuming SAME CHR!
##	--For probes that are not inside a feature, NA is returned instead.
##
##	f -- features;	chrom, chromStart, chromEnd
##	p -- probes;	chrom, chromStart, chromEnd
##
##	2009-05-23 Added to package groHMM from AffyTiling.
##	2009-02-06 Updated to allow assocation of features with arbitrary length.
##
########################################################################


associateWithInterval_foreachChrom <- function(i) {
	# Which KG?  prb?
	indxF   <- which(f[[1]] == C[i])
	indxPrb <- which(p[[1]] == C[i])

	if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
		# Type coersions.
		FeatureStart <- as.integer(f[[2]][indxF])
		FeatureEnd <- as.integer(f[[3]][indxF])
		PROBEStart <- as.integer(p[[2]][indxPrb])
		PROBELength <- as.integer(p[[3]][indxPrb] - p[[2]][indxPrb])

		# Set dimensions.
		dim(PROBEStart) <- c(NROW(PROBEStart), NCOL(PROBEStart))
		dim(FeatureStart) <- c(NROW(FeatureStart), NCOL(FeatureStart))
		dim(FeatureEnd) <- c(NROW(FeatureEnd), NCOL(FeatureEnd))
		dim(PROBELength) <- c(NROW(PROBELength), NCOL(PROBELength))

		Fprime <- .Call("AssociateRegionWithFeatures", FeatureStart, FeatureEnd, PROBEStart, PROBELength, PACKAGE = "groHMM")
	}
	
	return(Fprime)
}

#' associateWithInterval associates reads with the first feature that they fall inside of.
#'
#' @param features A GRanges object representing a set of genomic coordinates.  The meta-plot will be centered on the start position.
#' @param reads A GRanges object representing a set of mapped reads.
#' @return Returns index of the feature in which a each read is found.  Will be a vector of integers, the same size as the number of reads.  NA indicates that the reads does not fall inside of any feature.
#' @author Charles G. Danko and Minho Chae
associateWithInterval <- function(features, reads, ...) {
        f <- data.frame(chrom=as.character(seqnames(features)), start=as.integer(start(features)),
                                end=as.integer(end(features)), strand=as.character(strand(features)))
        p <- data.frame(chrom=as.character(seqnames(reads)), start=as.integer(start(reads)),
                                end=as.integer(end(reads)), strand=as.character(strand(reads))) 

	C <- as.character(unique(p[[1]]))
	mcp <- mclapply(c(1:NROW(C)), associateWithInterval_foreachChrom, ...)

	## Translate this into a return vector...
	F <- rep(NA, NROW(p))
	for(i in 1:NROW(C)) {
		indxPrb <- which(p[[1]] == C[i])
		F[indxPrb] <- as.character(f[[4]][indxF][as.vector(F[[i]])])
	}

    return(F)
}


