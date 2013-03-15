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
##	countReadsInInterval; also contains: countMappableReadsInInterval (below).
##	Date: 2009-05-07
##
##	This function takes two matricies -- both are essentially a BED file.
##		f  == features in a BED format, where columns represent: Chr, Start, End, Strand, ID.
##		p  == probes in a BED format, where columns represent: Chr, Start, End, Strand, ID.  
##   
##  Added 7/30/2010: Optinally, p can be a two-column data.frame with columns representing chrom, chromStart.  In this case strand is ignored.
##
##	Addumptions:
##	(1) FeatureStart and FeatureEnd should be vectors ONLY on the same 
##	chromosome as PROBEStart!!!
##	(2) Take care: only returns probes on the same strand as the feature, f.
##	to return all probes, force all strands to be the same.
##
##	TODO: 
##	(X) Write C function, which returns a vector representing the number of reads in each 
##	(2) ...
##
########################################################################

#' countReadsInInterval counts the number of reads that fall inside of each genomic feature.
#'
#' Only counts reads on the appropriate strand.
#'
#' @param features A GRanges object representing a set of genomic coordinates.  The meta-plot will be centered on the start position.
#' @param reads A GRanges object representing a set of mapped reads.
#' @return Returns a vector of counts, each representing the number of reads inside each genomic interval.
#' @author Charles G. Danko and Minho Chae
countReadsInInterval <- function(features, reads) {
	# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
	# Can't use GRanges sort b/c GRanges sort by natural order; order by (a) sequence level, (b)strand, 
	# (c)start, (d)width.

   	features <- features[order(as.character(seqnames(features)), start(features)),]
    reads <- reads[order(as.character(seqnames(reads)), start(reads)),]

	C <- sort(unique(as.character(seqnames(features))))
	F <- rep(0,NROW(features))
	for(i in 1:NROW(C)) {
		# Which KG?  prb?
		indxF   <- which(as.character(seqnames(features)) == C[i])
		indxPrb <- which(as.character(seqnames(reads)) == C[i])

		if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
			# Type coersions.
			FeatureStart 	<- start(features[indxF,])
			FeatureEnd 	<- end(features[indxF,])
			FeatureStr	<- as.character(strand(features[indxF,]))
			PROBEStart 	<- start(reads[indxPrb,])

		 	PROBEEnd        <- end(reads[indxPrb,])
                        PROBEStr        <- as.character(strand(reads[indxPrb,]))	

			#if(NCOL(p) > 2) { ## Assume that all four columns are present.
			#  PROBEEnd 	<- as.integer(p[indxPrb,3][Pord])
			#  PROBEStr	<- as.character(p[indxPrb,4][Pord])
			#}
			#else {
			# if(NCOL(p) == 2) { ## If probes are represented by only two columns, set PROBEEnd to 0-length probes.  Record 
			#   PROBEEnd <- PROBEStart
			#   PROBEStr <- rep("+", NROW(PROBEStart))
			#   FeatureStr <- rep("+", NROW(FeatureStr))
			# }
			#}

			# Set dimensions.
			dim(FeatureStart)	<- c(NROW(FeatureStart), NCOL(FeatureStart))
			dim(FeatureEnd) 	<- c(NROW(FeatureEnd), 	 NCOL(FeatureEnd))
			dim(FeatureStr)		<- c(NROW(FeatureStr), 	 NCOL(FeatureStr))
			dim(PROBEStart) 	<- c(NROW(PROBEStart), 	 NCOL(PROBEStart))
			dim(PROBEEnd) 		<- c(NROW(PROBEEnd), 	 NCOL(PROBEEnd))
			dim(PROBEStr)		<- c(NROW(PROBEStr), 	 NCOL(PROBEStr))

			Fprime <- .Call("CountReadsInFeatures", FeatureStart, FeatureEnd, FeatureStr, 
							PROBEStart, PROBEEnd, PROBEStr, PACKAGE = "groHMM")

			F[indxF] <- as.integer(Fprime)
		}
	}

	return(F)
}
 
########################################################################
##
##	countMappableReadsInInterval
##	Date: 2010-05-12
##
##	This function takes information from BED file to represent regions (as in CountReadsInInterval), and 
##	   a list structure representing unmappable positions.  Counts the number of mappable positions in 
##	   the interval.
##		f  == features in a BED format, where columns represent: Chr, Start, End, Strand, ID.
##		UnMAQ  == List structure of the un-mappable bases in the genome.
##
##	Addumptions:
##	(1) FeatureStart and FeatureEnd should be vectors ONLY on the same 
##	chromosome as PROBEStart!!!
##	(2) Take care: only returns probes on the same strand as the feature, f.
##	to return all probes, force all strands to be the same.
##
##	TODO: 
##	(1) ...
##	(2) ...
##
########################################################################

#' countMappableReadsInInterval counts the number of mappable reads in a set of genomic features.
#'
#' @param features A GRanges object representing a set of genomic coordinates.  The meta-plot will be centered on the start position.
#' @param UnMap List object representing the position of un-mappable reads.  Default: not used.
#' @param debug If set to TRUE, provides additional print options. Default: FALSE
#' @return Returns a vector of counts, each representing the number of reads inside each genomic interval.
#' @author Charles G. Danko and Minho Chae
countMappableReadsInInterval <- function(features, UnMap, debug=FALSE) {
	# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
    features <- features[order(as.character(seqnames(features)), start(features)),]

	C <- sort(unique(as.character(seqnames(features))))
	F <- rep(0,NROW(features))
	for(i in 1:NROW(C)) {
		indxF   <- which(as.character(seqnames(features)) == C[i])

		if(NROW(indxF) >0) {
			# Type coersions.
			FeatureStart 	<- start(features[indxF,])
			FeatureEnd 	<- end(features[indxF,])
			FeatureStr	<- as.character(strand(features[indxF,]))

			# Set dimensions.
			dim(FeatureStart)	<- c(NROW(FeatureStart), NCOL(FeatureStart))
			dim(FeatureEnd) 	<- c(NROW(FeatureEnd), 	 NCOL(FeatureEnd))
			dim(FeatureStr)		<- c(NROW(FeatureStr), 	 NCOL(FeatureStr))
			
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
				print(paste(C[i],": Counting unMAQable regions.",sep=""))
				print(paste("CHRSIZE:", CHRSIZE, "CHRSTART:", CHRSTART))
			}

			## Count unMAQable regions, and size of everything ... 
			nonmappable <- .Call("CountUnMAQableReads", FeatureStart, FeatureEnd, 
					UnMap[[2]], CHRSTART, CHRSIZE, PACKAGE = "groHMM")

			## Adjust size of gene body.
			Difference <- (FeatureEnd - FeatureStart) - nonmappable + 1 ## Otherwise, get -1 for some.

			if(debug) {
				print(head(nonmappable))
				print(as.integer(head(Difference)))
			}

			F[indxF] <- as.integer(Difference)
		}
	}

	return(F)
}
 
########################################################################
##
##  limitToXkb
##  Date: 2012-07-11
##
##  This function limits a genomic range to a samll region relative to the transcription site.
##
##
##  TODO:
##
########################################################################

#' limitToXkb truncates a set of genomic itnervals at a constant, maximum size.
#'
#' @param features A GRanges object representing a set of genomic coordinates.  The meta-plot will be centered on the start position.
#' @param offset Starts the interval from this position relative to the start of each genomic features.
#' @param size Specifies the size of the window.
#' @return Returns a new 'GRanges' object representing the size new size.
#' @author Minho Chae and Charles G. Danko
limitToXkb <- function(features, offset=1000, size=13000) {
	w <- width(features)

	# do nothing for w < offset 
	small  <- (offset < w) & (w < size)
	features[small,] <- flank(features[small,], -1*(w[small]-offset), start=FALSE)

	big  <- w > size 
	features[big,] <- resize(features[big,], width=size)

	bigPlus <- big & as.character(strand(features))=="+"
	if (any(bigPlus)) start(features[bigPlus,]) <- start(features[bigPlus,]) + offset 

	bigMinus <- big & as.character(strand(features))=="-"
	if (any(bigMinus)) end(features[bigMinus,]) <- end(features[bigMinus,]) - offset 

	return(features)
}


