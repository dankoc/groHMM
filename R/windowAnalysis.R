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
##	windowAnalysis
##	Date: 2009-07-08
##
##	Returns a vector of integers representing the counts of reads in a moving window.
##
##	Arguments:
##	p	-> data.frame of: CHR, START, END, STRAND.
##	strand	    -> Character; Strand of the probes to count, or "N" for ignore strand.
##	window_size	-> The size of the moving window.  Defaults to consecutive, non-inclusive windows.
##	step_size	-> The number of bp moved with each step.  Defaults to consecutive, non-inclusive windows.
##	chrom       -> Chromosome to search (NULL for all).
##	start       -> The start in genomic coordinates (NULL for all).
##	end	        -> The end in genomic coordinates (NULL for all).
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

windowAnalysis_foreachChrom <- function(i) {
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
		   PROBEStr <- as.character(rep(strand, NROW(PROBEStart))) ## no strand information takes both...
		 }
		 }

		if(limitPCRDups) {
		  print("WARNING: Using limitPCRDups assumes all probes are the same size!  Don't use for paired end data!!!!")
		  PROBElength <- (PROBEStart[1]-PROBEEnd[1])
		  PROBEStart <- as.integer(unique(PROBEStart))
		  PROBEEnd <- as.integer(PROBEStart+PROBElength)
		  PROBEStr <- as.character(rep(strand, NROW(PROBEStart)))
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
		Hprime <- .Call("WindowAnalysis", PROBEStart, PROBEEnd, PROBEStr, strand,
						wsize, ssize, start, endChrom, PACKAGE = "groHMM")

		#H[[chrom[i]]] <- as.integer(Hprime)
		return(as.integer(Hprime))
	}
	return(integer(0))
}


#' windowAnalysis Returns a vector of integers representing the counts of reads in a moving window.
#'
#' @param reads GenomicRanges object representing the position of reads mapping in the genome.
#' @param strand Takes values of "+", "-", or "N".  Computes Writes a wiggle on the speicified strand.  "N" denotes collapsing reads on both strands.  Default: "N".
#' @param window_size Size of the moving window. Either window_size or step_size must be specified.
#' @param step_size The number of bp moved with each step.
#' @param chrom Chromosome for which to return data.  Default: returns all avaliable data.
#' @param start The start position in genomic coordinates.  Default: returns all avaliable data.
#' @param end The end position in genomic coordinates.  Default: returns all avaliable data.
#' @param limitPCRDups Counts only one read mapping to each start site.  NOTE: If set to TRUE, assumes that all reads are the same length (don't use for paired-end data).  Default: FALSE.  
#' @param debug If set to TRUE, provides additional print options. Default: FALSE
#' @return List object, each element of which represents a chromosome.
#' @author Charles G. Danko and Minho Chae
windowAnalysis <- function(reads, strand="N", window_size=(step_size-1), step_size=(window_size+1), chrom=NULL, start=0, end=NULL, limitPCRDups=FALSE, debug=FALSE, ...) { 
	p <- data.frame(chrom=as.factor(as.character((seqnames(reads)))), start=as.integer(start(reads)),
                          end=as.integer(end(reads)), strand=as.factor(as.character(strand(reads))))

	wsize <- as.integer(window_size)
	ssize <- as.integer(step_size)
	strand   <- as.character(strand)
	start <- as.integer(start)
	
	if(limitPCRDups & strand != "N" & NCOL(p) > 2) {
	  p <- p[p[,4] == strand,]
	}
	if(is.null(chrom)) {
		chrom <- sort(as.character(unique(p[[1]])))
	}
	if(!is.null(end)) {
		endChrom <- as.integer(end)
	}

	H <- mclapply(c(1:NROW(chrom)), windowAnalysis_foreachChrom, ...)
	
	return(H)

}
 
