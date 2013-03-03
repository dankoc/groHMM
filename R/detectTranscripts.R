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


#######################################################################################################
##
##	detectTranscripts, detectTranscriptsUTS, detectTranscriptsV.
##	Date: 2009-07-17
##
##	Detects transcripts using a two-state HMM representing transcribed, and nontranscribed.
##
######################################################################################################


############################################################################
##
##	detectTranscriptsEM -- Runs full Baum-Welch EM to detect transcript
##	 positions and distribution paremeters in the same algorithm.
##
##	TODO: Make more general (for POLII ChIP-seq): 
##		strand="B", denotes both +/-; "N", denotes no information??
##
############################################################################

#' detectTranscripts detects transcripts de novo using a two-state hidden Markov model (HMM).
#'
#' Read counts can be specified as either a GRanges object (reads), or using a fixed-step wiggle-format passed in a list (Fp and Fm).  Either reads or BOTH Fp and Fm must be specified.
#'
#' @param features A GRanges object representing a set of genomic coordinates.  The meta-plot will be centered on the start position.
#' @param reads A GRanges object representing a set of mapped reads.
#' @param Fp Wiggle-formatted read counts on "+" strand. Optionally, Fp and Fm represent list() filled with a vector of counts for each chromosome.  Can detect transcripts starting from a fixed-step wiggle.
#' @param Fm Wiggle-formatted read counts on "-" strand. 
#' @param LtProbA Log probability of t... .  Default: -5. One of these is just an initialization, and the final value is set by EM.  The other is a holdout parameter.
#' @param LtProbB Log probability of t... .  Default: -5.
#' @param UTS Varience in read counts of the untranscribed sequence.  Default: 5.
#' @param size Log probability of t... .  Default: -5.
#' @param threshold Threshold change in total likelihood, below which EM exits. 
#' @return Returns a GRanges object representing the predicted genomic coordinates of transcripts on both the + and - strand.
#' @author Charles G. Danko and Minho Chae

## CGD: TODO: Test switch over to gamma, rather than dGamma?!
detectTranscripts <- function(reads=NULL, Fp=NULL, Fm=NULL, LtProbA=-5, LtProbB=-5, UTS=5, size=50, threshold=0.1, debug=TRUE) {

	stopifnot(!is.null(reads)|(!is.null(Fp) & !is.null(Fm)))

	## Setup/Window Analysis/Casting.
	epsilon <- 0.001
	
	if(is.null(Fp) & is.null(Fm)) { ## Allow equilavent form of Fp and Fm to be spcified in the function automatically.
	 Fp <- windowAnalysis(reads=reads, strand="+", ssize=size, debug=FALSE)
	 Fm <- windowAnalysis(reads=reads, strand="-", ssize=size, debug=FALSE)
	}
	
	nFp <- NROW(Fp)
	nFm <- NROW(Fm)
	CHRp <- as.character(names(Fp))
	CHRm <- as.character(names(Fm))

	size <- as.integer(size)
	ANS <- NULL

	## Set up initial HMM variables.
	HMM <- list()
	HMM$nstates <- as.integer(2)
	HMM$ePrDist <- c("dgamma", "dgamma") ## CGD: 3-3-13: Still legacy. Switch to integrating gamma between read and read+1

	HMM$iProb <- as.real(log(c(1.0,0.0)))
									## Non-transcribed,  transcribed.
	HMM$ePrVars <- as.list(data.frame(c(UTS, 1/UTS, -1), c(0.5, 10, -1)))
	HMM$tProb <- as.list(data.frame(c(log(1-exp(LtProbA)), LtProbA), c(LtProbB, log(1-exp(LtProbB))) ))

	## Cast counts to a real, and combine +/- strand into one list variable.  
	##  Treat like separate training sequences (they really are).
	F <- list()
	for(i in 1:nFp) F[[i]]     <- as.real(Fp[[i]]+1) ## CGD: 3-3-13: Still legacy.  Switch to integrating gamma between read and read+1
	for(i in 1:nFm) F[[i+nFp]] <- as.real(Fm[[i]]+1) ## CGD: 3-3-13: Still legacy.  Switch to integrating gamma between read and read+1

	## In case the above command copies, rather than points ... free unused memory.
	remove(Fp)
	remove(Fm)

	## Run EM algorithm.
	BWem <- .Call("RBaumWelchEM", HMM$nstates, F, as.integer(1),
				HMM$ePrDist, HMM$ePrVars, HMM$tProb, HMM$iProb, 
				as.real(threshold), c(TRUE, FALSE), c(FALSE, TRUE), as.integer(1), TRUE, PACKAGE="groHMM")
						# Update Transitions, Emissions.

	## Translate these into transcript positions.
	for(i in 1:NROW(CHRp)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i]]), as.real(0.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		# ANS <- rbind(ANS, data.frame(chrom =rep(CHRp[i], Nrep), chromStart =ans$Start, chromEnd =ans$End, 
		#	name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("+", Nrep)))
        	ANS <- rbind(ANS, data.frame(chrom =rep(CHRp[i], Nrep), start=ans$Start, end =ans$End, strand =rep("+", Nrep)))
	}

	for(i in 1:NROW(CHRm)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i+nFp]]), as.real(0.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		# ANS <- rbind(ANS, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)), chromStart =ans$Start, chromEnd =ans$End,
		# 	name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("-", Nrep)))
        	ANS <- rbind(ANS, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)), start=ans$Start, end=ans$End, strand =rep("-", Nrep)))
	}

	#BWem[[4]] <- ANS
	#names(BWem) <- c("EmisParams", "TransParams", "ViterbiStates", "Transcripts")

	BWem[[4]] <- GRanges(seqnames = Rle(ANS$chrom), ranges = IRanges(ANS$start, ANS$end-1), # ANS$end -1
                 		strand = Rle(strand(ANS$strand)), type=Rle("tx",NROW(ANS)), ID=paste(ANS$chrom, "_", ANS$start, ANS$strand, sep=""))
    	names(BWem) <- c("emisParams", "transParams", "viterbiStates", "transcripts")

	if(debug) {
		print(BWem[[1]])
		print(BWem[[2]])
	}

	return(BWem)
}
