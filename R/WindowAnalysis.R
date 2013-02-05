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
##	WindowAnalysis
##	Date: 2009-07-08
##
##	Returns a vector of integers representing the counts of reads in a moving window.
##
##	Arguments:
##	p	-> data.frame of: CHR, START, END, STRAND.
##	str	-> Character; Strand of the probes to count, or "N" for ignore strand.
##	wsize	-> The size of the moving window.  Defaults to consecutive, non-inclusive windows.
##	ssize	-> The number of bp moved with each step.  Defaults to consecutive, non-inclusive windows.
##	chrom   -> Chromosome to search (NULL for all).
##	start   -> The start in genomic coordinates (NULL for all).
##	end	-> The end in genomic coordinates (NULL for all).
##
##
##	Assumptions:
##	(1) 
##
##	TODO: 
##	(1) Write C function. 
##	(2) ...
##
########################################################################
WindowAnalysis <- function(p, str="N", wsize=(ssize-1), ssize=(wsize+1), chrom=NULL, start=0, end=NULL, limitPCRDups=FALSE, debug=FALSE) { 
	wsize <- as.integer(wsize)
	ssize <- as.integer(ssize)
	str   <- as.character(str)
	start <- as.integer(start)
	H <- NULL
	
	if(limitPCRDups & str != "N" & NCOL(p) > 2) {
	  p <- p[p[,4] == str,]
	}
	if(is.null(chrom)) {
		chrom <- sort(as.character(unique(p[[1]])))
	}
	if(!is.null(end)) {
		endChrom <- as.integer(end)
	}

	for(i in 1:NROW(chrom)) {
		if(debug) {
			print(chrom[i])
		}

		# Which KG?  prb?
		indxPrb <- which(as.character(p[[1]]) == chrom[i])

		if((NROW(indxPrb) >0)) {
			# Type coersions.
			PROBEStart 	<- as.integer(p[indxPrb,2])
			
			if(NCOL(p) > 2) { ## Assume that all four columns are present.
			  PROBEEnd 	<- as.integer(p[indxPrb,3])
			  PROBEStr	<- as.character(p[indxPrb,4])
			}
			else {
			 if(NCOL(p) == 2) { ## If probes are represented by only two columns, set PROBEEnd to 0-length probes.  Record 
			   PROBEEnd <- as.integer(PROBEStart)
			   PROBEStr <- as.character(rep(str, NROW(PROBEStart))) ## no strand information takes both...
			 }
			 }

			if(limitPCRDups) {
			  print("WARNING: Using limitPCRDups assumes all probes are the same size!!!!!")
			  PROBElength <- (PROBEStart[1]-PROBEEnd[1])
			  PROBEStart <- as.integer(unique(PROBEStart))
			  PROBEEnd <- as.integer(PROBEStart+PROBElength)
			  PROBEStr <- as.character(rep(str, NROW(PROBEStart)))
			}
			 
			# Set dimensions.
			dim(PROBEStart) 	<- c(NROW(PROBEStart), 	 NCOL(PROBEStart))
			dim(PROBEEnd) 		<- c(NROW(PROBEEnd), 	 NCOL(PROBEEnd))
			dim(PROBEStr)		<- c(NROW(PROBEStr), 	 NCOL(PROBEStr))

			if(is.null(end)) {
				endChrom <- as.integer(max(PROBEEnd))
			}

			if(debug) {
				print(paste(chrom[i],": Counting reads in specified region.",sep=""))
			}
			Hprime <- .Call("WindowAnalysis", PROBEStart, PROBEEnd, PROBEStr, str,
							wsize, ssize, start, endChrom, PACKAGE = "GROseq")

			H[[chrom[i]]] <- as.integer(Hprime)
		}
		if(debug) {
			print(paste(chrom[i],": Done!",sep=""))
		}
	}

	return(H)

}
 
