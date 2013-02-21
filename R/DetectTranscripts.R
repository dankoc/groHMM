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
##	detectTranscripts, detectTranscriptsEM, detectTranscriptsV.
##	Date: 2009-07-17
##
##	Detects transcripts using a two-state HMM representing transcribed, and nontranscribed.
##
######################################################################################################

## Sets up the HMM by translating the initList into vector over emission prob. for fitting paremeters.
setupInitialHMMVars <- function(initList, str, size, CHR, F, debug) {
#######	## Set up HMM vars...
	if(debug) print("Fitting assuming that initList represents transcripts.")
	nstates <- as.integer(2)
	ePrDist <- c("pois", "gamma")
	trans <- NULL
	unTrans <- NULL
	for(j in 1:NROW(CHR)) {

		## Increment to pseudocounts for the fit to gamma distribution (not 0!).
		F[[CHR[j]]] <- F[[CHR[j]]] + 1

		initListChr <- initList[which(initList[[1]] == CHR[j]),]
		initListChr[[2]] <- floor(initListChr[[2]]/size)
		initListChr[[3]] <- ceiling(initListChr[[3]]/size)

		if(NROW(initListChr) > 0) {
			## Translated regions.
			transNewIndx <- as.vector(unlist(lapply(c(1:NROW(initListChr)), function(X) { 
				c(initListChr[[2]][X]:initListChr[[3]][X]) })))
			transNewIndx <- unique(transNewIndx)
			transNewIndx <- transNewIndx[which(transNewIndx < NROW(F[[CHR[j]]]))]

			## Untranslated regions.
			untransNewIndx <- as.vector(unlist(lapply(c(2:NROW(initListChr)), function(X) { 
				c(initListChr[[3]][(X-1)]:initListChr[[2]][X]) })))
			untransNewIndx <- c(	c(1:initListChr[[2]][1]), 
						untransNewIndx, 
						c(initListChr[[3]][NROW(initListChr)]:NROW(F[[CHR[j]]])) )
			untransNewIndx <- unique(untransNewIndx)
			untransNewIndx <- trans[which(untransNewIndx < NROW(F[[CHR[j]]]))]

			## Catalogue emission probabilities
			trans <- c(trans,F[[CHR[j]]][transNewIndx])
			unTrans <- c(unTrans,F[[CHR[j]]][untransNewIndx])
		}
	}

#######	## Fit the emission probabilities.
	trans <- trans[which(is.finite(trans))]
	unTrans <- unTrans[which(is.finite(unTrans))]

	if(debug) {
		print(paste("mean:", mean(trans)))
		print(paste("var: ",  var(trans)))
	} 
	transcribedParams  <- RgammaMLE(trans)
	ePrVars <- as.list(data.frame(c(median(unTrans), -1, -1),#mean(unTrans), -1, -1),
				c(transcribedParams$shape, transcribedParams$scale, -1)))
	iProb <- as.real(log(c(1.0,0.0)))
	t_tr_un <- NROW(initList)/NROW(trans)
	t_un_tr <- NROW(initList)/NROW(unTrans)
	tProb <- as.list(data.frame(log(c(1-t_tr_un, t_tr_un)), log(c(t_un_tr, 1-t_un_tr))))


#######	## Write initial variables.  Very useful for debug!
	if(debug) {
		print(paste("NROW(Trans):", NROW(trans), "NROW(initList):", NROW(initList)))

		print("Emission Prob. Var.:")
		print(ePrVars)
		print("Initial Prob.:")
		print(iProb)
		print("Trans. Prob.:")
		print(tProb)
	}

######## Return paremeters as a list...
	returnList <- list()
	returnList$nstates <- nstates
	returnList$ePrDist <- ePrDist
	returnList$ePrVars <- ePrVars
	returnList$iProb <- iProb
	returnList$tProb <- tProb

	return(returnList)
}

######################################################################################################
##
##	detectTranscripts -- soon to be detectTranscriptsV
##
##	Arguments:
##	pgr	-> GRanges 
##	(deprecated: p	-> data.frame of: CHR, START, END, STRAND.)
##	initList-> BED file containing coordinates of training set (i.e. likely transcribed regions).
##	str	-> Character; Strand of the probes to count, or "N" for ignore strand.
##	size	-> The size of the moving window.  Defaults to consecutive, non-inclusive windows.
##
##
##	Addumptions:
##	(1) 
##
##	TODO: 
##	(1) For now ... just threshold raw reads ...
##	(2) Implement fitting to the background model, and transformation into -log10(P). 
##	(3) Use EM to detect transcripts ... 
##
######################################################################################################
detectTranscripts <- function(pgr, initList, strand="N", size=50, debug=TRUE) {
	## Setup/Window Analysis/Casting.
	F <- windowAnalysis(pgr=pgr, strand=strand, ssize=size, debug=FALSE)
	CHR <- as.character(names(F))
	size <- as.integer(size)
	ANS <- NULL

	HMM <- setupInitialHMMVars(initList, strand, size, CHR, F, debug)

	## Cast counts to a real.
	for(i in 1:NROW(F))
		F[[i]] <- as.real(F[[i]]+1)

#######	## Now run though and preduct using those variables...
	for(i in 1:NROW(CHR)) {
		## Run Viterbi algorithm.
		if(debug) {
			print(paste("Running Viterbi:", CHR[i],sep=" "))
		}
		viterbi <- .Call("Rviterbi", as.real(F[[CHR[i]]]), as.integer(NROW(F[[CHR[i]]])), as.integer(1), HMM$nstates, 
						HMM$ePrDist, HMM$ePrVars, HMM$tProb, HMM$iProb, PACKAGE="groHMM")

		## Use get transcript positions to convert state paths into BED files.
		if(debug) {
			print(paste("Getting Positions:", CHR[i],sep=" "))
		}
		ans <- .Call("getTranscriptPositions", as.real(viterbi), as.real(0.5), size, PACKAGE="groHMM")
		ANS <- rbind(ANS, data.frame(Chr =rep(CHR[i], NROW(ans$Start)), Start =ans$Start, End =ans$End))
	}

	return(ANS)
}

############################################################################
##
##	detectTranscriptsEM -- Runs full Baum-Welch EM to detect transcript
##	 positions and distribution paremeters in the same algorithm.
##
##	TODO: Make more general (for POLII ChIP-seq): 
##		strand="B", denotes both +/-; "N", denotes no information??
##
############################################################################
detectTranscriptsEM <- function(pgr=NULL, Fp=NULL, Fm=NULL, LtProbA=-5, LtProbB=-5, UTS=5, size=50, thresh=0.1, debug=TRUE) {

	stopifnot(!is.null(pgr)|(!is.null(Fp) & !is.null(Fm)))

	## Setup/Window Analysis/Casting.
	epsilon <- 0.001
	
	if(is.null(Fp) & is.null(Fm)) { ## Allow equilavent form of Fp and Fm to be spcified in the function automatically.
	 Fp <- windowAnalysis(pgr=pgr, strand="+", ssize=size, debug=FALSE)
	 Fm <- windowAnalysis(pgr=pgr, strand="-", ssize=size, debug=FALSE)
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
	HMM$ePrDist <- c("dgamma", "dgamma")

	HMM$iProb <- as.real(log(c(1.0,0.0)))
#	HMM$ePrVars <- as.list(data.frame(c(10, 1/10, -1), c(0.5, 10, -1)))
									## Non-transcribed,  transcribed.
	HMM$ePrVars <- as.list(data.frame(c(UTS, 1/UTS, -1), c(0.5, 10, -1)))
	HMM$tProb <- as.list(data.frame(c(log(1-exp(LtProbA)), LtProbA), c(LtProbB, log(1-exp(LtProbB))) ))

	## Cast counts to a real, and combine +/- strand into one list variable.  
	##  Treat like separate training sequences (they really are).
	F <- list()
	for(i in 1:nFp) F[[i]]     <- as.real(Fp[[i]]+1)
	for(i in 1:nFm) F[[i+nFp]] <- as.real(Fm[[i]]+1)

	## In case the above command copies, rather than points ... free unused memory.
	remove(Fp)
	remove(Fm)

	## Run EM algorithm.
	BWem <- .Call("RBaumWelchEM", HMM$nstates, F, as.integer(1),
				HMM$ePrDist, HMM$ePrVars, HMM$tProb, HMM$iProb, 
				as.real(thresh), c(TRUE, FALSE), c(FALSE, TRUE), as.integer(1), TRUE, PACKAGE="groHMM")
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


############################################################################
##
##	detectTranscriptsUTS -- Detects transcript positions using information
##                          from the plus and minus strand (groHMM data).
##                          The model is given in my notebook on 5/16/2011.
##
##  Arguments:
##         pgr       --> Sorted GRO-seq data in GRanges 
##         (deprecated: p         --> Sorted GRO-seq data (chrom, start, end, and strand).)
##         GENES     --> Sets the location of genes in teh dataset. (chrom, start, end, and strand).
##         ENH       --> Sets the location of putative enhancers in the dataset. (chrom and max).
##         UTS       --> Varience for the un-transcribed state.  Try setting this using EM?!
##
############################################################################
detectTranscriptsUTS <- function(pgr=NULL, Fp=NULL, Fm=NULL, GENES, ENH, UTS=5, size=50, thresh=0.1, debug=TRUE) {

    stopifnot(!is.null(pgr)|(!is.null(Fp) & !is.null(Fm)))

	## Setup/Window Analysis/Casting.
	epsilon <- 0.001
	
	if(is.null(Fp) & is.null(Fm)) { ## Allow equilavent form of Fp and Fm to be spcified in the function automatically.
	 Fp <- windowAnalysis(pgr=pgr, strand="+", ssize=size, debug=FALSE)
	 Fm <- windowAnalysis(pgr=pgr, strand="-", ssize=size, debug=FALSE)
	}
	
	nFp <- NROW(Fp)
	nFm <- NROW(Fm)
	CHRp <- as.character(names(Fp))
	CHRm <- as.character(names(Fm))
	stopifnot(sum(CHRp == CHRm)/NROW(CHRp) == 1) ## Sanity check -- require chromosomes be in the same order...
	
	## Cast counts to a real, and combine +/- strand into one list variable.
	##  Treat like separate training sequences (they really are).
	F <- list()
	for(i in 1:nFp) F[[i]]     <- as.real(Fp[[i]]+1)
	for(i in 1:nFm) F[[i+nFp]] <- as.real(Fm[[i]]+1)
	remove(Fp, Fm) ## In case the above command copies, rather than points ... free unused memory.

	size <- as.integer(size)
	ANS <- NULL

	## Set up initial HMM variables.
	HMM <- list()
	HMM$nstates <- as.integer(9) ## Nine states.
	HMM$iProb   <- as.real(log(c(1.0,rep(0.0, HMM$nstates-1))))

	##############################################################
	##
	## Emission prob.
	##
	HMM$nemis       <- as.integer(2) ## Two sets of emissions probabilities -- one per strand.
	HMM$updateEmis  <- rep(TRUE, HMM$nstates*HMM$nemis)
	HMM$ePrDist     <- rep("gamma", HMM$nstates*HMM$nemis)
	fitSize <- 250 #bp

	## States 2-3: Short, divergent, enhancer-like.  Use ~250bp upstream (- strand) and dowstream (+ strand).
	## Get position of enhancer. 
	sdelHigh <- NULL
	sdelLow  <- NULL
	for(c in unique(ENH[,1])) {
	 indx <- floor(ENH[ENH[,1]==c, 2]/size) ## Mid position divided by size should equate to the window index in the vector ...
	 sdelHigh <- c(sdelHigh, F[[which(CHRp==c)]][indx], F[[which(CHRm==c)+nFp]][indx])

	 for(count in 1:ceiling(fitSize/size)) {
	   sdelHigh <- c(sdelHigh, F[[which(CHRp==c)]][indx+count]) ## Plus strand...
	   sdelHigh <- c(sdelHigh,  F[[which(CHRm==c)+nFp]][indx-count]) ## Minus strand...

	   sdelLow  <- c(sdelLow,  F[[which(CHRp==c)]][indx-count]) ## Minus strand...
	   sdelLow  <- c(sdelLow,  F[[which(CHRm==c)+nFp]][indx+count]) ## Minus strand...
	 }
	}
	
	## Fit that vector to a gamma!!
	sdelHigh <- RgammaMLE(sdelHigh)
	sdelLow  <- RgammaMLE(sdelLow)
	
	## States 4-9: Genes ... consist of (A) divergent peak, (B) pause peak, and (C) gene body.
	## (A) divergent peak.  Use ~250bp (or fitSize) upstream of the TSS.
	## (B) pause peak.  Use ~250bp (or fitSize) downstream of the TSS.
	GENES <- GENES[(GENES[,3]-GENES[,2]) > (2*fitSize),]
	genesPauseHigh <- NULL
	genesPauseLow  <- NULL
	genesDivHigh   <- NULL
	genesDivLow    <- NULL
	genesGBHigh    <- NULL
	genesGBLow     <- NULL

	for(c in unique(GENES[,1])) {
	 indxp <- floor(GENES[GENES[,1]==c & GENES[,4]=="+", 2]/size) ## Get the index of the TSS...
	 indxm <- floor(GENES[GENES[,1]==c & GENES[,4]=="-", 3]/size)
	 genesPauseHigh <- c(genesPauseHigh, F[[which(CHRp==c)]][indxp], F[[which(CHRm==c)+nFp]][indxm])
	 genesDivHigh   <- c(genesDivHigh,   F[[which(CHRm==c)+nFp]][indxp], F[[which(CHRp==c)]][indxm])
	 
	 #for(count in 1:ceiling(fitSize/sz)) {  Minho: changed sz to size 2/11/13
	 for(count in 1:ceiling(fitSize/size)) { 
                                        ## Index for F represent strand of reads.  
                                                              ## indxp/m represent strand of the gene.
	   genesPauseHigh <- c(genesPauseHigh, F[[which(CHRp==c)]][indxp+count]) ## Genes on the plus strand...
	   genesPauseHigh <- c(genesPauseHigh,  F[[which(CHRm==c)+nFp]][indxm-count]) ## Minus strand... 
	   genesPauseLow  <- c(genesPauseLow,  F[[which(CHRm==c)+nFp]][indxp+count])
	   genesPauseLow  <- c(genesPauseLow,  F[[which(CHRp==c)]][indxm-count]) ## Minus strand...

	   genesDivHigh <- c(genesDivHigh, F[[which(CHRm==c)+nFp]][indxp-count]) 
	   genesDivHigh <- c(genesDivHigh,  F[[which(CHRp==c)]][indxm+count]) ## Genes on minus strand
	   genesDivLow  <- c(genesDivLow,  F[[which(CHRp==c)]][indxp-count])
	   genesDivLow  <- c(genesDivLow,  F[[which(CHRm==c)+nFp]][indxm+count]) ## Genes on minus strand
	 }

	 ## (C) gene body.
  	 GBindxP <- unique(unlist(lapply(c(1:sum(GENES[,1]==c&GENES[,4]=="+")),
		function(x) { floor((GENES[GENES[,1]==c&GENES[,4]=="+",2]+fitSize)/size)[x]:ceiling(GENES[GENES[,1]==c&GENES[,4]=="+",3]/size)[x]  })))
	 GBindxM <- unique(unlist(lapply(c(1:sum(GENES[,1]==c&GENES[,4]=="-")),
		function(x) { floor(GENES[GENES[,1]==c&GENES[,4]=="-",2]/size)[x]:ceiling((GENES[GENES[,1]==c&GENES[,4]=="-",3]-fitSize)/size)[x]  })))
	 genesGBHigh <- c(genesGBHigh, F[[which(CHRp==c)]][GBindxP])
	 genesGBHigh <- c(genesGBHigh, F[[which(CHRm==c)+nFp]][GBindxM])
	 genesGBLow  <- c(genesGBLow,  F[[which(CHRm==c)+nFp]][GBindxP])
	 genesGBLow  <- c(genesGBLow,  F[[which(CHRp==c)]][GBindxM])
	}
	
	## Fit that vector to a gamma!!
	genesPauseHigh<- RgammaMLE(genesPauseHigh)
	genesPauseLow <- RgammaMLE(genesPauseLow)
	genesDivHigh  <- RgammaMLE(genesDivHigh)
	genesDivLow   <- RgammaMLE(genesDivLow)
	genesGBHigh   <- RgammaMLE(genesGBHigh[!is.na(genesGBHigh)])
	genesGBLow    <- RgammaMLE(genesGBLow[!is.na(genesGBLow)])
	
										#################################
	                                    ## Plus strand emissions.
	HMM$ePrVars <- as.list(data.frame(c(UTS, 1/UTS, -1),                            ## Non-transcribed.

										## Short, divergent, enhancer-like.
										c(sdelLow$shape, sdelLow$scale, -1), c(sdelHigh$shape, sdelHigh$scale, -1),     

										## Plus strand transcript/gene: divergent, pause, body.
										c(genesDivLow$shape, genesDivLow$scale, -1), c(genesPauseHigh$shape, genesPauseHigh$scale, -1), c(genesGBHigh$shape, genesGBHigh$scale, -1),
										
										## Minus strand transcript: divergent, pause, body.
										c(genesDivHigh$shape, genesDivHigh$scale, -1), c(genesPauseLow$shape, genesPauseLow$scale, -1), c(genesGBLow$shape, genesGBLow$scale, -1),
										
										##################################
										## Minus strand emissions...
										c(UTS, 1/UTS, -1),							## Non-transcribed.

										## Short, divergent, enhancer-like.
										c(sdelHigh$shape, sdelHigh$scale, -1), c(sdelLow$shape, sdelLow$scale, -1), 	

										## Plus strand transcript: divergent, pause, body.
										c(genesDivHigh$shape, genesDivHigh$scale, -1), c(genesPauseLow$shape, genesPauseLow$scale, -1), c(genesGBLow$shape, genesGBLow$scale, -1),
										
										## Minus strand: divergent, pause, body.
										c(genesDivLow$shape, genesDivLow$scale, -1), c(genesPauseHigh$shape, genesPauseHigh$scale, -1), c(genesGBHigh$shape, genesGBHigh$scale, -1) ))

	##############################################################
	##
	## Transition prob.
	##																						## NOTE: FOR ALL STATES, SELF TRANSITIONS ALSO ALLOWED!
	HMM$updateTrans <- rep(TRUE, HMM$nstates)
	HMM$tProb <- as.list(log(data.frame(c(0.997, 0.001, 0, 0.001, 0, 0, 0, 0, 0.001),  			   ## State 1 (NonTranscribed) --> 2, 4, 9
									c(0, 0.95, 0.05, 0, 0, 0, 0, 0, 0),  						   ## State 2 (short, divergent) ---> 3
									c(0.05, 0, 0.95, 0, 0, 0, 0, 0, 0),                            ## State 3 (short, pause) ---> 1
									c(0, 0, 0, 0.95, 0.05, 0, 0, 0, 0),                            ## State 4 (trans. plus, div) ---> 5
									c(0, 0, 0, 0, 0.95, 0.05, 0, 0, 0),                            ## State 5 (trans. plus, pause) ---> 6
									c(0.001, 0, 0, 0, 0, 0.999, 0, 0, 0),                          ## State 6 (trans. plus, body) ---> 1
									c(0.05, 0, 0, 0, 0, 0, 0.95, 0, 0),                            ## State 7 (trans. minus, div) ---> 1
									c(0, 0, 0, 0, 0, 0, 0.05, 0.95, 0),                            ## State 8 (trans. minus, pause) ---> 7
									c(0, 0, 0, 0, 0, 0, 0, 0.001, 0.999))*50/size))                ## State 9 (trans. minus, body) ---> 8

	##############################################################
	##
	## Run EM algorithm.
	##
	BWem <- .Call("RBaumWelchEM", HMM$nstates, F, HMM$nemis,
				HMM$ePrDist, HMM$ePrVars, HMM$tProb, HMM$iProb, 
				as.real(thresh), HMM$updateTrans, HMM$updateEmis, as.integer(1), TRUE, PACKAGE="groHMM")
						# Update Transitions, Emissions.

	## Translate these into transcript positions.
	for(i in 1:NROW(CHRp)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i]]), as.real(0.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		ANS <- rbind(ANS, data.frame(chrom =rep(CHRp[i], Nrep), chromStart =ans$Start, chromEnd =ans$End, 
			name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("+", Nrep)))
	}

	for(i in 1:NROW(CHRm)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i+nFp]]), as.real(0.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		ANS <- rbind(ANS, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)), chromStart =ans$Start, chromEnd =ans$End,
			name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("-", Nrep)))
	}

	BWem[[4]] <- ANS
	names(BWem) <- c("EmisParams", "TransParams", "ViterbiStates", "Transcripts")

	if(debug) {
		print(BWem[[1]])
		print(BWem[[2]])
	}

	return(BWem)
}

############################################################################
##
##	DetectDecreaseIncrease -- Identifies systematic increase in the GRO-seq data.
##
##  Arguments:
##         pgr         --> Sorted GRO-seq data in GRanges
##         (deprecated:	p         --> Sorted GRO-seq data (chrom, start, end, and strand).
##         GENES     --> Sets the location of genes in teh dataset. (chrom, start, end, and strand).
##         size      --> Size of the sliding window (bp).
##         localRegion   --> The half-life of the upstream/downstream window used to form the background information (bp).
##
############################################################################

detectDecreaseIncrease <- function(pgr=NULL, Fp=NULL, Fm=NULL, size=50, localRegion=50000, debug=TRUE) {
	stopifnot(!is.null(pgr)|(!is.null(Fp) & !is.null(Fm)))
	if(is.null(Fp) & is.null(Fm)) { ## Allow equilavent form of Fp and Fm to be spcified in the function automatically.
	 Fp <- windowAnalysis(pgr=pgr, strand="+", ssize=size, debug=FALSE)
	 Fm <- windowAnalysis(pgr=pgr, strand="-", ssize=size, debug=FALSE)
	}

	nFp <- NROW(Fp)
	nFm <- NROW(Fm)
	CHRp <- as.character(names(Fp))
	CHRm <- as.character(names(Fm))
	stopifnot(sum(CHRp == CHRm)/NROW(CHRp) == 1) ## Sanity check -- require chromosomes be in the same order...

	F <- list()
	for(i in 1:nFp) F[[i]]     <- as.real(Fp[[i]])
	for(i in 1:nFm) F[[i+nFp]] <- as.real(Fm[[i]])
	remove(Fp, Fm) ## In case the above command copies, rather than points ... free unused memory.

	###############################
	## Get the weights... Take ~20 half-lives.
	Lambda <- size/localRegion
	wgt <- pexp(c(0:round(20/Lambda)), rate=Lambda, lower.tail=FALSE)
	nBlocksAnalyzed <- NROW(wgt)
	
	peaksTrack <- list()
	for(chr in c(1:NROW(F))) {
		if(debug) print(paste("Analyzing chrom",chr,"of",NROW(F)))
		nWind <- NROW(F[[chr]]) ## Size of the chromsome (in windows).
		peaksTrack[[chr]] <- c(0,unlist(lapply(c(2:(nWind-1)), function(x) {
			indxL<- c((x-nBlocksAnalyzed):(x-1))
			indxR<- c((x+1):(x+nBlocksAnalyzed))
			
			if(chr > nFp) {
				indxUpstream <- indxR
				indxDownstream <- indxL
			} ## If - strand, upstream is right.
			else {
				indxUpstream <- rev(indxL)
				indxDownstream <- rev(indxR)
			} ## else (+ strand), upstream is left.  Reverse it to be in the same order as the weights.
			indx <- indx[(indx>0)&(indx<=NROW(F[[chr]]))] ## NOTE: Could result in 0... and/ or LOTs of na's.  
			wgtIteration <- wgt[c(1:NROW(indx))] ## Because on the 5' or 3' end of the chromsome, indx may not be as long as teh list of weights.
			
			mu<- weighted.mean(x= F[[chr]][indx], w= wgtIteration, na.rm=TRUE) #mean(F[[chr]][indx], na.rm=TRUE) ## OR: use weighted mean...
			if(is.na(mu) | mu < muGlobal) mu <- muGlobal
			va<- weighted.var(x= F[[chr]][indx], w= wgtIteration, na.rm=TRUE) ##var(F[[chr]][indx], na.rm=TRUE)
			
			if(!is.na(va) & va>mu) { ## if overdispersed -- use NB.
			 sz<-(mu^2)/abs(va - mu)
			 return(-1*pnbinom(q=F[[chr]][x], mu=mu, size=sz, log.p=TRUE, lower.tail=FALSE))
			}
			else { ## else, use Poisson.
			  return(-1*ppois(q=F[[chr]][x], lambda= mu, log.p=TRUE, lower.tail=FALSE))
			}
		})),0)
		peaksTrack[[chr]][is.na(peaksTrack[[chr]])] <- 0
		if(debug) print(head(peaksTrack[[chr]]))
	}

	return(peaksTrack)
}

############################################################################
##
##	detectPeaks -- Identifies peaks in the GRO-seq data.
##
##  Arguments:
##         pg        --> Sorted GRO-seq data in GRanges 
##         (deprecated: p         --> Sorted GRO-seq data (chrom, start, end, and strand).)
##         size      --> Size of the sliding window (bp).
##         localRegion   --> The half-life of the upstream/downstream window used to form the background information (bp).
##
##  Returns: list comprised of:
##			[[1]]    --> Peaks
##			[[2]]    --> Decrease
##          [[3]]    --> Increase
##
############################################################################

detectPeaks <- function(pgr=NULL, Fp=NULL, Fm=NULL, size=50, localRegion=10000, debug=TRUE) {

	stopifnot(!is.null(pgr)|(!is.null(Fp) & !is.null(Fm)))
	if(is.null(Fp) & is.null(Fm)) { ## Allow equilavent form of Fp and Fm to be spcified in the function automatically.
	 Fp <- windowAnalysis(pgr=pgr, strand="+", ssize=size, debug=FALSE)
	 Fm <- windowAnalysis(pgr=pgr, strand="-", ssize=size, debug=FALSE)
	}

	nFp <- NROW(Fp)
	nFm <- NROW(Fm)
	CHRp <- as.character(names(Fp))
	CHRm <- as.character(names(Fm))
	stopifnot(sum(CHRp == CHRm)/NROW(CHRp) == 1) ## Sanity check -- require chromosomes be in the same order...

	F <- list()
	for(i in 1:nFp) F[[i]]     <- as.real(Fp[[i]])
	for(i in 1:nFm) F[[i+nFp]] <- as.real(Fm[[i]])
	remove(Fp, Fm) ## In case the above command copies, rather than points ... free unused memory.

	###############################
	## Evaluate peak-caller.
	upstreamSize <- localRegion/ size ## Size of the radius to search(in windows).

	dv <- real()
	for(chr in c(1:NROW(F))) dv <- c(dv, F[[chr]])
	muGlobal <- mean(dv, na.rm=TRUE)
	vaGlobal <- var(dv, na.rm=TRUE)
	szGlobal <- (muGlobal^2)/(var(dv, na.rm=TRUE) - muGlobal)
	if(debug) print(paste("global mean=",muGlobal,"size=",szGlobal,"varience=",vaGlobal))
	remove(dv)

	###############################
	## Get the weights... Take ~20 half-lives.
	Lambda <- size/localRegion
	wgt <- pexp(c(0:round(20/Lambda)), rate=Lambda, lower.tail=FALSE)
	nBlocksAnalyzed <- NROW(wgt)
	
	allTracks <- list()
	peaksTrack <- list()
	tssTrack <- list()
	changeTrack <- list()
	changeWilcox <- list()
	for(chr in c(1:NROW(F))) {
		if(debug) print(paste("Analyzing chrom",chr,"of",NROW(F)))
		nWind <- NROW(F[[chr]]) ## Size of the chromsome (in windows).
		dC <- unlist(lapply(c(2:(nWind-1)), function(x) {
			indxL<- c((x-nBlocksAnalyzed):(x-1))
			indxR<- c((x+1):(x+nBlocksAnalyzed))
			
			if(chr > nFp) {
				indx <- indxR
				indxDns <- rev(indxL)
			} ## If - strand, upstream is right.
			else {
				indx <- rev(indxL)
				indxDns <- indxR
			} ## else (+ strand), upstream is left.  Reverse it to be in the same order as the weights.
			indx <- indx[(indx>0)&(indx<=NROW(F[[chr]]))] ## NOTE: Could result in 0... and/ or LOTs of na's.  
			indxDns <- indxDns[(indxDns>0)&(indxDns<=NROW(F[[chr]]))] ## NOTE: Could result in 0... and/ or LOTs of na's.  
			wgtIteration <- wgt[c(1:NROW(indx))] ## Because on the 5' or 3' end of the chromsome, indx may not be as long as teh list of weights.
			wgtIterationDns <- wgt[c(1:NROW(indxDns))] ## Because on the 5' or 3' end of the chromsome, indx may not be as long as teh list of weights.

			mu<- weighted.mean(x= F[[chr]][indx], w= wgtIteration, na.rm=TRUE) #mean(F[[chr]][indx], na.rm=TRUE) ## OR: use weighted mean...
			muDns<- weighted.mean(x= F[[chr]][indxDns], w= wgtIterationDns, na.rm=TRUE) #mean(F[[chr]][indx], na.rm=TRUE) ## OR: use weighted mean...
			muBoth<- weighted.mean(x= F[[chr]][c(indx,indxDns)], w= c(wgtIteration,wgtIterationDns), na.rm=TRUE) #mean(F[[chr]][indx], na.rm=TRUE) ## OR: use weighted mean...

			if(is.na(mu) | mu < muGlobal) mu <- muGlobal
			if(is.na(muDns) | muDns < muGlobal) muDns <- muGlobal
			if(is.na(muBoth) | muBoth < muGlobal) muBoth <- muGlobal

			va<- weighted.var(x= F[[chr]][indx], w= wgtIteration, na.rm=TRUE) ##var(F[[chr]][indx], na.rm=TRUE)
			vaDns<- weighted.var(x= F[[chr]][indxDns], w= wgtIterationDns, na.rm=TRUE) ##var(F[[chr]][indx], na.rm=TRUE)
			vaBoth<- weighted.var(x= F[[chr]][c(indx,indxDns)], w= c(wgtIteration,wgtIterationDns), na.rm=TRUE) ##var(F[[chr]][indx], na.rm=TRUE)
		
			## Null model is that downstream counts fit w/ the upstream distribution paramters (mu, sz).  Alternative model is that it fits better with downstram params (muDns, szDns).
			sz<-(mu^2)/abs(va - mu)
			szDns<-(muDns^2)/abs(vaDns - muDns)
			szBoth<-(muBoth^2)/abs(vaBoth - muBoth)
			
			if(!is.na(vaBoth) & vaBoth>muBoth) { ## if overdispersed -- use NB.
			 pVpeak <- -1*pnbinom(q=F[[chr]][x], mu=muBoth, size=szBoth, log.p=TRUE, lower.tail=FALSE)
			 tss <- -1*pnbinom(q=F[[chr]][x], mu=mu, size=sz, log.p=TRUE, lower.tail=FALSE)
 			 lrtDl <- -2*sum(wgtIterationDns*dnbinom(x=F[[chr]][indxDns], mu=mu, size=sz, log=TRUE)) +2*sum(wgtIterationDns*dnbinom(x=F[[chr]][indxDns], mu=muDns, size=szDns, log=TRUE))
  			 lrtDr <- -2*sum(wgtIteration*dnbinom(x=F[[chr]][indx], mu=muDns, size=szDns, log=TRUE)) +2*sum(wgtIteration*dnbinom(x=F[[chr]][indx], mu=mu, size=sz, log=TRUE))
			}
			else { ## else, use Poisson.
			  pVpeak <- -1*ppois(q=F[[chr]][x], lambda= muBoth, log.p=TRUE, lower.tail=FALSE)
			  tss <- -1*ppois(q=F[[chr]][x], lambda= mu, log.p=TRUE, lower.tail=FALSE)
			  lrtDl <- -2*sum(wgtIterationDns*dpois(x=F[[chr]][indxDns], lambda=mu, log=TRUE)) +2*sum(wgtIterationDns*dpois(x=F[[chr]][indxDns], lambda=muDns, log=TRUE))
			  lrtDr <- -2*sum(wgtIteration*dpois(x=F[[chr]][indx], lambda=muDns, log=TRUE)) +2*sum(wgtIteration*dpois(x=F[[chr]][indx], lambda=mu, log=TRUE))
			}
			hcUps <- c(1:min(NROW(indx), round(localRegion/size)))
			hcDns <- c(1:min(NROW(indxDns), round(localRegion/size)))
			difMWUpsLess <- 	-log(wilcox.test(x= F[[chr]][indx[hcUps]], y= F[[chr]][indxDns[hcDns]], alternative="greater")$p.value)
			difMWDnsLess <- 	-log(wilcox.test(x= F[[chr]][indx[hcUps]], y= F[[chr]][indxDns[hcDns]], alternative="less")$p.value)
			lrtD <- (which.max(c(lrtDr,0,lrtDl))-2)*max(lrtDr, lrtDl)
			difMW <- (which.max(c(difMWUpsLess,0,difMWDnsLess))-2)*max(difMWUpsLess, difMWDnsLess)
			
		    return(c(pVpeak,tss,lrtD,difMW))
		}))

		peaksTrack[[chr]] <- c(0,dC[seq(1,NROW(dC),4)],0)
		peaksTrack[[chr]][is.na(peaksTrack[[chr]])] <- 0
		
		tssTrack[[chr]] <- c(0,dC[seq(2,NROW(dC),4)],0)
		tssTrack[[chr]][is.na(tssTrack[[chr]])] <- 0
		
		changeTrack[[chr]] <- c(0,dC[seq(3,NROW(dC),4)],0)
		changeTrack[[chr]][is.na(changeTrack[[chr]])] <- 0
		
		changeWilcox[[chr]] <- c(0,dC[seq(4,NROW(dC),4)],0)
		changeWilcox[[chr]][is.na(changeWilcox[[chr]])] <- 0
		
		if(debug) { 
			print(head(peaksTrack[[chr]])) 
			print(head(changeTrack[[chr]]))
		}
	}
	allTracks[[1]] <- peaksTrack
	allTracks[[2]] <- tssTrack
	allTracks[[3]] <- changeTrack
	allTracks[[4]] <- changeWilcox
	names(allTracks) <- c("peaksTrack", "tssTrack", "ChangeNB", "ChangeWilcox")
	return(allTracks)
}

############################################################################
##
##	detectTranscriptsPeaks -- Detects transcript positions and enhancers/ TSS
##                            using information on which windows contain peaks 
##                            in the groHMM data. The model is given in my 
##                            notebook on 5/17/2011 and 5/19/2011.
##
##  Arguments:
##         pgr       --> Sorted GRO-seq data in GRanges 
##         (deprecatd: p         --> Sorted GRO-seq data (chrom, start, end, and strand).)
##         GENES     --> Sets the location of genes in teh dataset. (chrom, start, end, and strand).
##         ENH       --> Sets the location of putative enhancers in the dataset. (chrom and max).
##         UTS       --> Varience for the un-transcribed state.  Try setting this using EM?!
##
############################################################################
detectTranscriptsPeaks <- function(pgr=NULL, Fp=NULL, Fm=NULL, GENES, ENH, UTS=5, size=50, thresh=0.1, debug=TRUE) {

    stopifnot(!is.null(pgr)|(!is.null(Fp) & !is.null(Fm)))

	## Setup/Window Analysis/Casting.
	epsilon <- 0.001
	
	if(is.null(Fp) & is.null(Fm)) { ## Allow equilavent form of Fp and Fm to be spcified in the function automatically.
	 Fp <- windowAnalysis(pgr=pgr, strand="+", ssize=size, debug=FALSE)
	 Fm <- windowAnalysis(pgr=pgr, strand="-", ssize=size, debug=FALSE)
	}
	
	nFp <- NROW(Fp)
	nFm <- NROW(Fm)
	CHRp <- as.character(names(Fp))
	CHRm <- as.character(names(Fm))
	stopifnot(sum(CHRp == CHRm)/NROW(CHRp) == 1) ## Sanity check -- require chromosomes be in the same order...

	##############################################################
	##
	## Evaluate peak-caller.
	##
	upstreamSize <- 10000/size ## 10kb?!
	unlist(lapply( , function(x) {
		dpois()
	}))

hist(Fp[[1]][floor(GENES[14,2]/size):ceiling(GENES[14,3]/size)])
	
	## Cast counts to a real, and combine +/- strand into one list variable.
	##  Treat like separate training sequences (they really are).
	F <- list()
	for(i in 1:nFp) F[[i]]     <- as.real(Fp[[i]]+1)
	for(i in 1:nFm) F[[i+nFp]] <- as.real(Fm[[i]]+1)
	remove(Fp, Fm) ## In case the above command copies, rather than points ... free unused memory.

	size <- as.integer(size)
	ANS <- NULL

	## Set up initial HMM variables.
	HMM <- list()
	HMM$nstates <- as.integer(9) ## Nine states.
	HMM$iProb   <- as.real(log(c(1.0,rep(0.0, HMM$nstates-1))))

	##############################################################
	##
	## Emission prob.
	##
	HMM$nemis       <- as.integer(2) ## Two sets of emissions probabilities -- one per strand.
	HMM$updateEmis  <- rep(TRUE, HMM$nstates*HMM$nemis)
	HMM$ePrDist     <- rep("gamma", HMM$nstates*HMM$nemis)
	fitSize <- 500 #bp

	## States 2-3: Short, divergent, enhancer-like.  Use ~250bp upstream (- strand) and dowstream (+ strand).
	## Get position of enhancer. 
	sdelHigh <- NULL
	sdelLow  <- NULL
	for(c in unique(ENH[,1])) {
	 indx <- floor(ENH[ENH[,1]==c, 2]/size) ## Mid position divided by size should equate to the window index in the vector ...
	 sdelHigh <- c(sdelHigh, F[[which(CHRp==c)]][indx], F[[which(CHRm==c)+nFp]][indx])

	 for(count in 1:ceiling(fitSize/size)) {
	   sdelHigh <- c(sdelHigh, F[[which(CHRp==c)]][indx+count]) ## Plus strand...
	   sdelHigh <- c(sdelHigh,  F[[which(CHRm==c)+nFp]][indx-count]) ## Minus strand...

	   sdelLow  <- c(sdelLow,  F[[which(CHRp==c)]][indx-count]) ## Minus strand...
	   sdelLow  <- c(sdelLow,  F[[which(CHRm==c)+nFp]][indx+count]) ## Minus strand...
	 }
	}
	
	## Fit that vector to a gamma!!
	sdelHigh <- RgammaMLE(sdelHigh)
	sdelLow  <- RgammaMLE(sdelLow)
	
	## States 4-9: Genes ... consist of (A) divergent peak, (B) pause peak, and (C) gene body.
	## (A) divergent peak.  Use ~250bp (or fitSize) upstream of the TSS.
	## (B) pause peak.  Use ~250bp (or fitSize) downstream of the TSS.
	GENES <- GENES[(GENES[,3]-GENES[,2]) > (2*fitSize),]
	genesPauseHigh <- NULL
	genesPauseLow  <- NULL
	genesDivHigh   <- NULL
	genesDivLow    <- NULL
	genesGBHigh    <- NULL
	genesGBLow     <- NULL

	for(c in unique(GENES[,1])) {
	 indxp <- floor(GENES[GENES[,1]==c & GENES[,4]=="+", 2]/size) ## Get the index of the TSS...
	 indxm <- floor(GENES[GENES[,1]==c & GENES[,4]=="-", 3]/size)
	 genesPauseHigh <- c(genesPauseHigh, F[[which(CHRp==c)]][indxp], F[[which(CHRm==c)+nFp]][indxm])
	 genesDivHigh   <- c(genesDivHigh,   F[[which(CHRm==c)+nFp]][indxp], F[[which(CHRp==c)]][indxm])
	 
	 for(count in 1:ceiling(fitSize/size)) {
                                        ## Index for F represent strand of reads.  
                                                              ## indxp/m represent strand of the gene.
	   genesPauseHigh <- c(genesPauseHigh, F[[which(CHRp==c)]][indxp+count]) ## Genes on the plus strand...
	   genesPauseHigh <- c(genesPauseHigh,  F[[which(CHRm==c)+nFp]][indxm-count]) ## Minus strand... 
	   genesPauseLow  <- c(genesPauseLow,  F[[which(CHRm==c)+nFp]][indxp+count])
	   genesPauseLow  <- c(genesPauseLow,  F[[which(CHRp==c)]][indxm-count]) ## Minus strand...

	   genesDivHigh <- c(genesDivHigh, F[[which(CHRm==c)+nFp]][indxp-count]) 
	   genesDivHigh <- c(genesDivHigh,  F[[which(CHRp==c)]][indxm+count]) ## Genes on minus strand
	   genesDivLow  <- c(genesDivLow,  F[[which(CHRp==c)]][indxp-count])
	   genesDivLow  <- c(genesDivLow,  F[[which(CHRm==c)+nFp]][indxm+count]) ## Genes on minus strand
	 }

	 ## (C) gene body.
  	 GBindxP <- unique(unlist(lapply(c(1:sum(GENES[,1]==c&GENES[,4]=="+")),
		function(x) { floor((GENES[GENES[,1]==c&GENES[,4]=="+",2]+fitSize)/size)[x]:ceiling(GENES[GENES[,1]==c&GENES[,4]=="+",3]/size)[x]  })))
	 GBindxM <- unique(unlist(lapply(c(1:sum(GENES[,1]==c&GENES[,4]=="-")),
		function(x) { floor(GENES[GENES[,1]==c&GENES[,4]=="-",2]/size)[x]:ceiling((GENES[GENES[,1]==c&GENES[,4]=="-",3]-fitSize)/size)[x]  })))
	 genesGBHigh <- c(genesGBHigh, F[[which(CHRp==c)]][GBindxP])
	 genesGBHigh <- c(genesGBHigh, F[[which(CHRm==c)+nFp]][GBindxM])
	 genesGBLow  <- c(genesGBLow,  F[[which(CHRm==c)+nFp]][GBindxP])
	 genesGBLow  <- c(genesGBLow,  F[[which(CHRp==c)]][GBindxM])
	}
	
	## Fit that vector to a gamma!!
	genesPauseHigh<- RgammaMLE(genesPauseHigh)
	genesPauseLow <- RgammaMLE(genesPauseLow)
	genesDivHigh  <- RgammaMLE(genesDivHigh)
	genesDivLow   <- RgammaMLE(genesDivLow)
	genesGBHigh   <- RgammaMLE(genesGBHigh[!is.na(genesGBHigh)])
	genesGBLow    <- RgammaMLE(genesGBLow[!is.na(genesGBLow)])
	
										#################################
	                                    ## Plus strand emissions.
	HMM$ePrVars <- as.list(data.frame(c(UTS, 1/UTS, -1),                            ## Non-transcribed.

										## Short, divergent, enhancer-like.
										c(sdelLow$shape, sdelLow$scale, -1), c(sdelHigh$shape, sdelHigh$scale, -1),     

										## Plus strand transcript/gene: divergent, pause, body.
										c(genesDivLow$shape, genesDivLow$scale, -1), c(genesPauseHigh$shape, genesPauseHigh$scale, -1), c(genesGBHigh$shape, genesGBHigh$scale, -1),
										
										## Minus strand transcript: divergent, pause, body.
										c(genesDivHigh$shape, genesDivHigh$scale, -1), c(genesPauseLow$shape, genesPauseLow$scale, -1), c(genesGBLow$shape, genesGBLow$scale, -1),
										
										##################################
										## Minus strand emissions...
										c(UTS, 1/UTS, -1),							## Non-transcribed.

										## Short, divergent, enhancer-like.
										c(sdelHigh$shape, sdelHigh$scale, -1), c(sdelLow$shape, sdelLow$scale, -1), 	

										## Plus strand transcript: divergent, pause, body.
										c(genesDivHigh$shape, genesDivHigh$scale, -1), c(genesPauseLow$shape, genesPauseLow$scale, -1), c(genesGBLow$shape, genesGBLow$scale, -1),
										
										## Minus strand: divergent, pause, body.
										c(genesDivLow$shape, genesDivLow$scale, -1), c(genesPauseHigh$shape, genesPauseHigh$scale, -1), c(genesGBHigh$shape, genesGBHigh$scale, -1) ))

	##############################################################
	##
	## Transition prob.
	##																						## NOTE: FOR ALL STATES, SELF TRANSITIONS ALSO ALLOWED!
	HMM$updateTrans <- rep(TRUE, HMM$nstates)
	HMM$tProb <- as.list(log(data.frame(c(0.997, 0.001, 0, 0.001, 0, 0, 0, 0, 0.001),  			   ## State 1 (NonTranscribed) --> 2, 4, 9
									c(0, 0.95, 0.05, 0, 0, 0, 0, 0, 0),  						   ## State 2 (short, divergent) ---> 3
									c(0.05, 0, 0.95, 0, 0, 0, 0, 0, 0),                            ## State 3 (short, pause) ---> 1
									c(0, 0, 0, 0.95, 0.05, 0, 0, 0, 0),                            ## State 4 (trans. plus, div) ---> 5
									c(0, 0, 0, 0, 0.95, 0.05, 0, 0, 0),                            ## State 5 (trans. plus, pause) ---> 6
									c(0.001, 0, 0, 0, 0, 0.999, 0, 0, 0),                          ## State 6 (trans. plus, body) ---> 1
									c(0.05, 0, 0, 0, 0, 0, 0.95, 0, 0),                            ## State 7 (trans. minus, div) ---> 1
									c(0, 0, 0, 0, 0, 0, 0.05, 0.95, 0),                            ## State 8 (trans. minus, pause) ---> 7
									c(0, 0, 0, 0, 0, 0, 0, 0.001, 0.999))*50/size))                ## State 9 (trans. minus, body) ---> 8

	##############################################################
	##
	## Run EM algorithm.
	##
	BWem <- .Call("RBaumWelchEM", HMM$nstates, F, HMM$nemis,
				HMM$ePrDist, HMM$ePrVars, HMM$tProb, HMM$iProb, 
				as.real(thresh), HMM$updateTrans, HMM$updateEmis, TRUE, PACKAGE="groHMM")
						# Update Transitions, Emissions.

	## Translate these into transcript positions.
	for(i in 1:NROW(CHRp)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i]]), as.real(0.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		ANS <- rbind(ANS, data.frame(chrom =rep(CHRp[i], Nrep), chromStart =ans$Start, chromEnd =ans$End, 
			name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("+", Nrep)))
	}

	for(i in 1:NROW(CHRm)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i+nFp]]), as.real(0.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		ANS <- rbind(ANS, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)), chromStart =ans$Start, chromEnd =ans$End,
			name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("-", Nrep)))
	}

	BWem[[4]] <- ANS
	names(BWem) <- c("EmisParams", "TransParams", "ViterbiStates", "Transcripts")

	if(debug) {
		print(BWem[[1]])
		print(BWem[[2]])
	}

	return(BWem)
}

############################################################################
##
##	DetectTranscriptsTSSPeakCalls -- Detects transcripts using the TSS Peak Calls as an additional covariate.
##
##	TODO: Make more general (for POLII ChIP-seq): 
##		strand="B", denotes both +/-; "N", denotes no information??
##
############################################################################
detectTranscriptsTSSPeakCalls <- function(pgr=NULL, Fp=NULL, Fm=NULL, TSSEmis, LtProbA=-5, LtProbB=-50, UTS=5, size=100, thresh=0.1, debug=TRUE) {

    stopifnot(!is.null(pgr)|(!is.null(Fp) & !is.null(Fm)))

	## Setup/Window Analysis/Casting.
	epsilon <- 0.001
	
	if(is.null(Fp) & is.null(Fm)) { ## Allow equilavent form of Fp and Fm to be spcified in the function automatically.
	 Fp <- windowAnalysis(pgr=pgr, strand="+", ssize=size, debug=FALSE)
	 Fm <- windowAnalysis(pgr=pgr, strand="-", ssize=size, debug=FALSE)
	}
	
	nFp <- NROW(Fp)
	nFm <- NROW(Fm)
	CHRp <- as.character(names(Fp))
	CHRm <- as.character(names(Fm))

	size <- as.integer(size)
	ANS <- NULL
	ANS2 <- NULL

	## Set up initial HMM variables.
	HMM <- list()
	HMM$nstates <- as.integer(3)
	HMM$nemis	<- as.integer(3)
	HMM$ePrDist <- c("gamma", "dgamma", "dgamma", "gamma", "dgamma", "dgamma", "gamma", "dgamma", "dgamma")
	HMM$updateEmis <- c(FALSE, TRUE, TRUE)
	HMM$updateTrans<- c(TRUE, TRUE, FALSE)

	HMM$iProb <- as.real(log(c(1.0,0.0,0.0)))
#	HMM$ePrVars <- as.list(data.frame(c(10, 1/10, -1), c(0.5, 10, -1)))
									## Non-transcribed,  peak, transcribed.
	HMM$ePrVars <- as.list(data.frame(c(UTS, 1/UTS, -1), c(0.5, 10, -1), c(0.5, 10, -1), ## Reads.
									  c(UTS, 1/UTS, -1), c(1, 30, -1), c(0.5, 5, -1),    ## Peaks1.
									  c(UTS, 1/UTS, -1), c(1, 30, -1), c(0.5, 5, -1)))   ## Peaks2.
	HMM$tProb <- as.list(data.frame(c(log(1-exp(LtProbA)), LtProbA, 0), 
									c(2*LtProbA, LtProbA, log(1-3*exp(LtProbA))),
									c(LtProbB, 2*LtProbB, log(1-3*exp(LtProbB)))))

	## Cast counts to a real, and combine +/- strand into one list variable.  
	##  Treat like separate training sequences (they really are).
	F <- list()
	for(i in seq(0,nFp-1,1)) {
		F[[i*HMM$nemis+1]]     <- as.real(Fp[[i+1]]+0.5)
		F[[i*HMM$nemis+2]]   <- as.real(TSSEmis[[i+1]]) ###### DOUBLE CHECK INDEXING!
		F[[i*HMM$nemis+3]]   <- as.real(TSSEmis[[i+2]]) ###### DOUBLE CHECK INDEXING!
	}
	for(i in seq(0,nFm-1,1)) {
		F[[i*HMM$nemis+nFp*HMM$nemis+1]] <- as.real(Fm[[i+1]]+0.5)
		F[[i*HMM$nemis+nFp*HMM$nemis+2]] <- as.real(TSSEmis[[i+2]]) ###### DOUBLE CHECK INDEXING!
		F[[i*HMM$nemis+nFp*HMM$nemis+3]] <- as.real(TSSEmis[[i+1]]) ###### DOUBLE CHECK INDEXING!
	}
	for(i in c(1:NROW(F))) {
		print(NROW(F[[i]]))
	}

	## In case the above command copies, rather than points ... free unused memory.
	remove(Fp)
	remove(Fm)

	## Run EM algorithm.
	BWem <- .Call("RBaumWelchEM", HMM$nstates, F, HMM$nemis,
				HMM$ePrDist, HMM$ePrVars, HMM$tProb, HMM$iProb, 
				as.real(thresh), HMM$updateTrans, HMM$updateEmis, as.integer(1), TRUE, PACKAGE="groHMM")
						# Update Transitions, Emissions.

	## Translate these into transcript positions.
	for(i in 1:NROW(CHRp)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i]]), as.real(0.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		ANS <- rbind(ANS, data.frame(chrom =rep(CHRp[i], Nrep), chromStart =ans$Start, chromEnd =ans$End, 
			name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("+", Nrep)))
	}

	for(i in 1:NROW(CHRm)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i+nFp]]), as.real(0.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		ANS <- rbind(ANS, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)), chromStart =ans$Start, chromEnd =ans$End,
			name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("-", Nrep)))
	}
	
	for(i in 1:NROW(CHRp)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i]]), as.real(1.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		ANS2 <- rbind(ANS2, data.frame(chrom =rep(CHRp[i], Nrep), chromStart =ans$Start, chromEnd =ans$End, 
			name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("+", Nrep)))
	}
	for(i in 1:NROW(CHRm)) {
		ans <- .Call("getTranscriptPositions", as.real(BWem[[3]][[i+nFp]]), as.real(1.5), size, PACKAGE="groHMM")
		Nrep <- NROW(ans$Start)
		ANS2 <- rbind(ANS2, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)), chromStart =ans$Start, chromEnd =ans$End,
			name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("-", Nrep)))
	}


	BWem[[4]] <- ANS
	BWem[[5]] <- ANS2
	names(BWem) <- c("EmisParams", "TransParams", "ViterbiStates", "Transcripts_peakAndBody", "Bodies")

	if(debug) {
		print(BWem[[1]])
		print(BWem[[2]])
	}

	return(BWem)
}


