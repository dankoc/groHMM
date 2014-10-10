###########################################################################
##
##   Copyright 2013, 2014 Charles Danko and Minho Chae.
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


#' detectTranscripts detects transcripts de novo using a two-state hidden 
#' Markov model (HMM).
#'
#' Read counts can be specified as either a GRanges object (reads), or using a 
#' fixed-step wiggle-format passed in a list (Fp and Fm).  
#' Either reads or BOTH Fp and Fm must be specified.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors set the option 'mc.cores'.
#'
#'  Reference: Hah N, Danko CG, Core L, Waterfall JJ, Siepel A, Lis JT, 
#' Kraus WL. A rapid, extensive, and transient transcriptional response to 
#' estrogen signaling in breast cancer cells. Cell. 2011 May 13;145(4):622-34. 
#' doi: 10.1016/j.cell.2011.03.042. 
#'
#' @param reads A GRanges object representing a set of mapped reads.
#' @param Fp Wiggle-formatted read counts on "+" strand. Optionally, Fp and Fm 
#' represent list() filled with a vector of counts for each chromosome.  
#' Can detect transcripts starting from a fixed-step wiggle.
#' @param Fm Wiggle-formatted read counts on "-" strand. 
#' @param LtProbA Log probability of t... .  Default: -5. One of these is just 
#' an initialization, and the final value is set by EM.  The other is a holdout
#' parameter.
#' @param LtProbB Log probability of t... .  Default: -200.
#' @param UTS Varience in read counts of the untranscribed sequence.  
#' Default: 5.
#' @param size Log probability of t... .  Default: -5.
#' @param threshold Threshold change in total likelihood, below which EM exits.
#' @param debug If set to TRUE, provides additional print options. 
#' Default: FALSE
#' @param ... Extra argument passed to mclapply
#' @return Returns a list of emisParams, trnasParams, viterbiStates, and 
#' transcripts.  The transcript element is a GRanges object representing the 
#' predicted genomic coordinates of transcripts on both the + and - strand.
#' @author Charles G. Danko and Minho Chae
#' @examples
#' S0mR1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam",
#'                package="groHMM")), "GRanges")
#' ## Not run:
#' # hmmResult <- detectTranscripts(S0mR1, LtProbB=-200, UTS=5, threshold=1)
#' # txHMM <- hmmResult$transcripts
## CGD: TODO: Test switch over to gamma, rather than dGamma?!

detectTranscripts <- function(reads=NULL, Fp=NULL, Fm=NULL, LtProbA=-5, 
    LtProbB=-200, UTS=5, size=50, threshold=0.1, debug=TRUE, ...) {

    stopifnot(!is.null(reads)|(!is.null(Fp) & !is.null(Fm)))

    ## Setup/Window Analysis/Casting.
    epsilon <- 0.001

    ## Allow equilavent form of Fp and Fm to be spcified in the function 
    ## automatically.
    if(is.null(Fp) & is.null(Fm)) { 
     Fp <- windowAnalysis(reads=reads, strand="+", windowSize=size, ...)
     Fm <- windowAnalysis(reads=reads, strand="-", windowSize=size, ...)
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
    ## CGD: 3-3-13: Still legacy. Switch to integrating gamma between read 
    ## and read+1

    HMM$iProb <- as.double(log(c(1.0,0.0)))
                                    ## Non-transcribed,  transcribed.
    HMM$ePrVars <- as.list(data.frame(c(UTS, 1/UTS, -1), c(0.5, 10, -1)))
    HMM$tProb <- as.list(data.frame(c(log(1-exp(LtProbA)), LtProbA), 
        c(LtProbB, log(1-exp(LtProbB))) ))

    ## Cast counts to a real, and combine +/- strand into one list variable.  
    ##  Treat like separate training sequences (they really are).
    FT <- list()    # MHC; 7/2/2014, bioconductor complains F thinking as False
    for(i in 1:nFp) FT[[i]]<- as.double(Fp[[i]]+1) 
    ## CGD: 3-3-13: Still legacy.  Switch to integrating gamma between read and 
    ## read+1

    for(i in 1:nFm) FT[[i+nFp]] <- as.double(Fm[[i]]+1) 
    ## CGD: 3-3-13: Still legacy.  Switch to integrating gamma between read and
    ## read+1

    ## In case the above command copies, rather than points ... 
    ## free unused memory.
    remove(Fp)
    remove(Fm)

    ## Run EM algorithm.
    BWem <- .Call("RBaumWelchEM", HMM$nstates, FT, as.integer(1),
                HMM$ePrDist, HMM$ePrVars, HMM$tProb, HMM$iProb, 
                as.double(threshold), c(TRUE, FALSE), c(FALSE, TRUE), 
                as.integer(1), TRUE, PACKAGE="groHMM")
               # Update Transitions, Emissions.

    ## Translate these into transcript positions.
    for(i in seq_along(CHRp)) {
        ans <- .Call("getTranscriptPositions", as.double(BWem[[3]][[i]]), 
                    as.double(0.5), size, PACKAGE="groHMM")
        Nrep <- NROW(ans$Start)
        # ANS <- rbind(ANS, data.frame(chrom =rep(CHRp[i], Nrep), 
        #           chromStart =ans$Start, chromEnd =ans$End, 
        #   name =rep("N", Nrep), score =rep("1", Nrep), strand =rep("+", 
        #           Nrep)))
            ANS <- rbind(ANS, data.frame(chrom =rep(CHRp[i], Nrep), 
                start=ans$Start, end =ans$End, strand =rep("+", Nrep)))
    }

    for(i in seq_along(CHRm)) {
        ans <- .Call("getTranscriptPositions", as.double(BWem[[3]][[i+nFp]]), 
                    as.double(0.5), size, PACKAGE="groHMM")
        Nrep <- NROW(ans$Start)
        # ANS <- rbind(ANS, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)), 
        #           chromStart =ans$Start, chromEnd =ans$End, 
        #           name =rep("N", Nrep), score =rep("1", Nrep), 
        #           strand =rep("-", Nrep)))
            ANS <- rbind(ANS, data.frame(chrom =rep(CHRm[i], NROW(ans$Start)), 
                    start=ans$Start, end=ans$End, strand =rep("-", Nrep)))
    }

    #BWem[[4]] <- ANS
    #names(BWem) <- c("EmisParams", "TransParams", "ViterbiStates", 
    #                   "Transcripts")

    BWem[[4]] <- GRanges(seqnames = Rle(ANS$chrom), 
                    ranges = IRanges(ANS$start, ANS$end-1), 
                    strand = Rle(strand(ANS$strand)), type=Rle("tx",NROW(ANS)), 
                    ID=paste(ANS$chrom, "_", ANS$start, ANS$strand, sep=""))
        names(BWem) <- c("emisParams", "transParams", "viterbiStates", 
                            "transcripts")

    if(debug) {
        print(BWem[[1]])
        print(BWem[[2]])
    }

    return(BWem)
}
