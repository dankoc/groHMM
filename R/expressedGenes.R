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
##	expressedGenes
##	Date: 2010-02-23
##
##	This identifes genes that are expressed in a given cell, based on short read data.
##		f  == genes/annotations; columns represent: Chr, Start, End, Strand, ID.
##		p  == short reads; columns represent: Chr, Start, End, Strand, ID.
##		UnMAQ == unmappable regions in the genome of interest.
##
##	Function defines expression as in Core, Waterfall, Lis; Science, Dec. 2008.
##
##
########################################################################

#' Function identifies expressed features using the methods introduced in Core, Waterfall, Lis; Science, Dec. 2008.
#'
#' @param features A GRanges object representing a set of genomic coordinates.  The meta-plot will be centered on the start position.
#' @param reads A GRanges object representing a set of mapped reads.
#' @param genomeSize The size of the target genome.  Default: 3e9, or roughly the size of the human genome.
#' @param Lambda Measurement of assay noise.  Default: # reads/ genome size (tends to be too high for GRO-seq data).
#' @param UnMap List object representing the position of un-mappable reads.  Default: not used.
#' @param debug If set to true, returns the number of positions.  Default: FALSE.
#' @return A data.frame representing the expression p.values for features of interest.
#' @author Charles G. Danko and Minho Chae
expressedGenes <- function(features, reads, genomeSize=3e9, Lambda= NULL, UnMap=NULL, debug=FALSE) {
	C <- sort(as.character(unique(features[[1]])))
	ANSgeneid <- rep("char", NROW(features))
	ANSpvalue <- rep(0,NROW(features))
	ANScounts <- rep(0,NROW(features))
	ANSunmapp <- rep(0,NROW(features))
	ANSgsize  <- rep(0,NROW(features))
	if(is.null(Lambda))	Lambda <- NROW(reads)/genomeSize
	for(i in 1:NROW(C)) {
		if(debug) {
			print(paste("Doing chromosome", C[i]))
		}
	
		# Which KG?  prb?
		indxF   <- which(as.character(features[[1]]) == C[i])
		indxPrb <- which(as.character(reads[[1]]) == C[i])

		if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
			# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
			Ford <- order(features[indxF,2])
			Pord <- order(reads[indxPrb,2])

			# Type coersions.
			FeatureStart 	<- as.integer(features[indxF,2][Ford])
			FeatureEnd 	<- as.integer(features[indxF,3][Ford])
			FeatureStr	<- as.character(features[indxF,4][Ford])
			PROBEStart 	<- as.integer(reads[indxPrb,2][Pord])
			PROBEEnd 	<- as.integer(reads[indxPrb,3][Pord])
			PROBEStr	<- as.character(reads[indxPrb,4][Pord])

			# Set dimensions.
			dim(FeatureStart)	<- c(NROW(FeatureStart), NCOL(FeatureStart))
			dim(FeatureEnd) 	<- c(NROW(FeatureEnd), 	 NCOL(FeatureEnd))
			dim(FeatureStr)		<- c(NROW(FeatureStr), 	 NCOL(FeatureStr))
			dim(PROBEStart) 	<- c(NROW(PROBEStart), 	 NCOL(PROBEStart))
			dim(PROBEEnd) 		<- c(NROW(PROBEEnd), 	 NCOL(PROBEEnd))
			dim(PROBEStr)		<- c(NROW(PROBEStr), 	 NCOL(PROBEStr))

			
			NUMReads <- .Call("CountReadsInFeatures", FeatureStart, FeatureEnd, FeatureStr, 
							PROBEStart, PROBEEnd, PROBEStr, PACKAGE = "groHMM")


			## Calculate UN-MAQable regions...
			if(!is.null(UnMap)) {

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
				MappablePositions <- (FeatureEnd - FeatureStart) - nonmappable + 1

				if(debug) {
					print(summary(nonmappable))
					print(summary(MappablePositions))
				}
			}
			else {
				nonmappable <- 1
			}
			
			MappablePositions <- (FeatureEnd - FeatureStart) - nonmappable + 1

			## Calculate poisson prob. of each.
			ANSgeneid[indxF][Ford] <- as.character(features[indxF, 5][Ford])
			ANSpvalue[indxF][Ford] <- ppois(NUMReads, (Lambda*MappablePositions), lower.tail=FALSE)
			ANScounts[indxF][Ford] <- NUMReads
			ANSunmapp[indxF][Ford] <- nonmappable
			ANSgsize[indxF][Ford]  <- (FeatureEnd-FeatureStart)
		}
	}

	return(data.frame(ID= ANSgeneid, pval= ANSpvalue, readCounts= ANScounts, 
					nonMappablePositions= ANSunmapp, size= ANSgsize))
}
 
