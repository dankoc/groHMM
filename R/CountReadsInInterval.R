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
##	CountReadsInInterval; also contains: CountMappableReadsInInterval (below).
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
CountReadsInInterval <- function(f, p) {
 
	C <- sort(as.character(unique(f[[1]])))
	F <- rep(0,NROW(f))
	for(i in 1:NROW(C)) {
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
			if(NCOL(p) > 2) { ## Assume that all four columns are present.
			  PROBEEnd 	<- as.integer(p[indxPrb,3][Pord])
			  PROBEStr	<- as.character(p[indxPrb,4][Pord])
			}
			else {
			 if(NCOL(p) == 2) { ## If probes are represented by only two columns, set PROBEEnd to 0-length probes.  Record 
			   PROBEEnd <- PROBEStart
			   PROBEStr <- rep("+", NROW(PROBEStart))
			   FeatureStr <- rep("+", NROW(FeatureStr))
			 }
			}

			# Set dimensions.
			dim(FeatureStart)	<- c(NROW(FeatureStart), NCOL(FeatureStart))
			dim(FeatureEnd) 	<- c(NROW(FeatureEnd), 	 NCOL(FeatureEnd))
			dim(FeatureStr)		<- c(NROW(FeatureStr), 	 NCOL(FeatureStr))
			dim(PROBEStart) 	<- c(NROW(PROBEStart), 	 NCOL(PROBEStart))
			dim(PROBEEnd) 		<- c(NROW(PROBEEnd), 	 NCOL(PROBEEnd))
			dim(PROBEStr)		<- c(NROW(PROBEStr), 	 NCOL(PROBEStr))

			Fprime <- .Call("CountReadsInFeatures", FeatureStart, FeatureEnd, FeatureStr, 
							PROBEStart, PROBEEnd, PROBEStr, PACKAGE = "GROseq")

			F[indxF][Ford] <- as.integer(Fprime)
		}
	}

	return(F)
}
 
########################################################################
##
##	CountMappableReadsInInterval
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
CountMappableReadsInInterval <- function(f, UnMAQ, debug=FALSE) {

	C <- sort(as.character(unique(f[[1]])))
	F <- rep(0,NROW(f))
	for(i in 1:NROW(C)) {
		indxF   <- which(as.character(f[[1]]) == C[i])

		if(NROW(indxF) >0) {
			# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
			Ford <- order(f[indxF,2])

			# Type coersions.
			FeatureStart 	<- as.integer(f[indxF,2][Ford])
			FeatureEnd 	<- as.integer(f[indxF,3][Ford])
			FeatureStr	<- as.character(f[indxF,4][Ford])

			# Set dimensions.
			dim(FeatureStart)	<- c(NROW(FeatureStart), NCOL(FeatureStart))
			dim(FeatureEnd) 	<- c(NROW(FeatureEnd), 	 NCOL(FeatureEnd))
			dim(FeatureStr)		<- c(NROW(FeatureStr), 	 NCOL(FeatureStr))
			
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
					UnMAQ[[2]], CHRSTART, CHRSIZE, PACKAGE = "GROseq")

			## Adjust size of gene body.
			Difference <- (FeatureEnd - FeatureStart) - nonmappable + 1 ## Otherwise, get -1 for some.

			if(debug) {
				print(head(nonmappable))
				print(as.integer(head(Difference)))
			}

			F[indxF][Ford] <- as.integer(Difference)
		}
	}

	return(F)
}
 
