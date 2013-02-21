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
expressedGenes <- function(f, p, genomeSize=3000000000, Lambda= NULL, UnMAQ=NULL, debug=FALSE) {
	C <- sort(as.character(unique(f[[1]])))
	ANSgeneid <- rep("char", NROW(f))
	ANSpvalue <- rep(0,NROW(f))
	ANScounts <- rep(0,NROW(f))
	ANSunmapp <- rep(0,NROW(f))
	ANSgsize  <- rep(0,NROW(f))
	if(is.null(Lambda))	Lambda <- NROW(p)/genomeSize
	for(i in 1:NROW(C)) {
		if(debug) {
			print(paste("Doing chromosome", C[i]))
		}
	
		# Which KG?  prb?
		indxF   <- which(as.character(f[[1]]) == C[i])
		indxPrb <- which(as.character(p[[1]]) == C[i])

		if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
			# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
			Ford <- order(f[indxF,2])
			Pord <- order(p[indxPrb,2])

			# Type coersions.
			FeatureStart 	<- as.integer(f[indxF,2][Ford])
			FeatureEnd 	<- as.integer(f[indxF,3][Ford])
			FeatureStr	<- as.character(f[indxF,4][Ford])
			PROBEStart 	<- as.integer(p[indxPrb,2][Pord])
			PROBEEnd 	<- as.integer(p[indxPrb,3][Pord])
			PROBEStr	<- as.character(p[indxPrb,4][Pord])

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
			if(!is.null(UnMAQ)) {

				## Count start index.
				chr_indx <- which(UnMAQ[[1]][[1]] == C[i])
				CHRSIZE <- as.integer(UnMAQ[[1]][[2]][chr_indx])
				CHRSTART <- as.integer(0)
				if(chr_indx > 1) {  ## Running on 1:0 gives c(1, 0)
					CHRSTART <- as.integer( 
						sum(UnMAQ[[1]][[2]][
							c(1:(chr_indx-1))
						]) +1)
				}

				if(debug) {
					print(paste(C[i],": Counting unMAQable regions.",sep=""))
					print(paste("CHRSIZE:", CHRSIZE, "CHRSTART:", CHRSTART))
				}

				## Count unMAQable regions, and size of everything ... 
				nonmappable <- .Call("CountUnMAQableReads", FeatureStart, FeatureEnd, 
						UnMAQ[[2]], CHRSTART, CHRSIZE, PACKAGE = "groHMM")

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
			ANSgeneid[indxF][Ford] <- as.character(f[indxF, 5][Ford])
			ANSpvalue[indxF][Ford] <- ppois(NUMReads, (Lambda*MappablePositions), lower.tail=F)
			ANScounts[indxF][Ford] <- NUMReads
			ANSunmapp[indxF][Ford] <- nonmappable
			ANSgsize[indxF][Ford]  <- (FeatureEnd-FeatureStart)
		}
	}

	return(data.frame(ID= ANSgeneid, pval= ANSpvalue, readCounts= ANScounts, 
					nonMappablePositions= ANSunmapp, size= ANSgsize))
}
 
