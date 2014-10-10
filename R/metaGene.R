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


#' Returns a histogram of the number of reads in each section of a moving 
#' window centered on a certain feature.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors, set the option 'mc.cores'.
#'
#' @param features A GRanges object representing a set of genomic coordinates.
#' The meta-plot will be centered on the transcription start site (TSS)
#' @param reads A GRanges object representing a set of mapped reads.  
#' Instead of 'reads', 'plusCVG' and 'minusCVG' can be used  Default: NULL
#' @param plusCVG A RangesList object for reads with '+' strand. 
#' @param minusCVG A RangesList object for reads with '-' strand. 
#' @param size The size of the moving window.
#' @param up Distance upstream of each features to align and histogram. 
#' Default: 10 kb.
#' @param down Distance downstream of each features to align and histogram. 
#' If NULL, same as up. Default: NULL.
#' @param ... Extra argument passed to mclapply
#' @return Returns a integer-Rle representing the 'typical' signal 
#' centered on a point of interest.
#' @author Charles G. Danko and Minho Chae
#' @examples
#' features <- GRanges("chr7", IRanges(1000, 1000), strand="+")
#' reads <- GRanges("chr7", IRanges(start=c(1000:1004, 1100), 
#'  width=rep(1, 6)), strand="+")
#' mg <- metaGene(features, reads, size=4, up=10)
metaGene <- function(features, reads=NULL, plusCVG=NULL, minusCVG=NULL, 
    size=100L, up=10000L, down=NULL, ...) {
    seqlevels(features) <- seqlevelsInUse(features)
    # Check 'reads'
    if (is.null(reads)) {
        if (is.null(plusCVG) || is.null(minusCVG))
            stop("Either 'reads' or 'plusCVG' and 'minusCVG' must be used")
    } else {
        seqlevels(reads) <- seqlevelsInUse(reads)
        plusCVG <- coverage(reads[strand(reads)=="+",])
        minusCVG <- coverage(reads[strand(reads)=="-",])
    }
    if (is.null(down)) down <- up

    featureList <- split(features, seqnames(features))

    H <- mclapply(seqlevels(features), metaGene_foreachChrom, 
        featureList=featureList, plusCVG=plusCVG, minusCVG=minusCVG, 
        size=size, up=up, down=down, ...)
    M <- sapply(seq_len(length(H)), function(x) as.integer(H[[x]]))

    return(Rle(apply(M, 1, sum)))
}


metaGene_foreachChrom <- function(chrom, featureList, plusCVG, minusCVG, 
    size, up, down) {
        f <- featureList[[chrom]]

        pCVG <- plusCVG[[chrom]]
        mCVG <- minusCVG[[chrom]]
        offset <- floor(size/2L)

        pro <- promoters(f, upstream=up+offset, downstream=(down+offset-1L))

        M <- sapply(1:length(pro), function(x) {
                if (as.character(strand(pro)[x]) == "+")
                    as.integer(runsum(pCVG[start(pro)[x]:end(pro)[x]], 
                        k=size))
                else
                    as.integer(rev(runsum(mCVG[start(pro)[x]:end(pro)[x]], 
                        k=size)))
            })
        return(Rle(apply(M, 1, sum)))
}


#' Runs metagene analysis for sense and antisense direction.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors, set the option 'mc.cores'.
#'
#' @param features GRanges A GRanges object representing a set of genomic 
#' coordinates, i.e., set of genes.
#' @param reads GRanges of reads.
#' @param anchorType Either 'TSS' or 'TTS'.  Metagene will be centered on the 
#' transcription start site(TSS) or transcription termination site(TTS).  
#' Default: TSS.
#' @param size Numeric.  The size of the moving window. Default: 100L
#' @param normCounts Numeric.  Normalization vector such as average reads.  
#' Default: 1L
#' @param up Numeric. Distance upstream of each feature to align and histogram.
#' Default: 1 kb
#' @param down Numeric. Distance downstream of each feature to align and 
#' histogram.  If NULL, down is same as up. Default: NULL
#' @param sampling Logical.  If TRUE, subsampling of Metagene is used.  
#' Default: FALSE
#' @param nSampling Numeric. Number of subsampling.  Default: 1000L
#' @param samplingRatio Numeric. Ratio of sampling for features.  Default: 0.1
#' @param ... Extra argument passed to mclapply.
#' @return A list of integer-Rle for sense and antisene.
#' @author Minho Chae
#' @examples
#' features <- GRanges("chr7", IRanges(start=1000:1001, width=rep(1,2)), 
#'  strand=c("+", "-"))
#' reads <- GRanges("chr7", IRanges(start=c(1000:1003, 1100:1101), 
#'  width=rep(1, 6)), strand=rep(c("+","-"), 3))
#' ## Not run:
#' # mg <- runMetaGene(features, reads, size=4, up=10)
runMetaGene <- function(features, reads, anchorType="TSS", size=100L, 
    normCounts=1L, up=10000L, down=NULL, sampling=FALSE, nSampling=1000L,
    samplingRatio=0.1, ...) {
    # Check 'anchorType'
    if (!anchorType %in% c("TSS", "TTS")) {
        stop("'anchorType' must be either 'TSS' or 'TTS'")
    }

    if (anchorType == "TSS") {
        f <- resize(features, width=1L, fix="start")
    } else if (anchorType == "TTS") {
        f <- resize(features, width=1L, fix="end")
    }

    if (is.null(down)) down <- up

    fRev <- f
    strand(fRev) <- rev(strand(f))

    plusCVG <- coverage(reads[strand(reads)=="+",])
    minusCVG <- coverage(reads[strand(reads)=="-",])
    
    message("sense ... ", appendLF=FALSE) 
    if (sampling) { 
        sense <- samplingMetaGene(features=f, plusCVG=plusCVG, 
            minusCVG=minusCVG, size=size, up=up, down=down, 
            nSampling=nSampling, samplingRatio=samplingRatio, ...) 
    } else {
        sense <- metaGene(features=f, plusCVG=plusCVG, minusCVG=minusCVG, 
                    size=size, up=up, down=down, ...)
        sense <- sense/length(features)
    }
    message("OK")

    message("antisense ... ", appendLF=FALSE)
    if (sampling) {
        antisense <- samplingMetaGene(features=fRev, plusCVG=plusCVG, 
            minusCVG=minusCVG, size=size, up=up, down=down,
            nSampling=nSampling, samplingRatio=samplingRatio, ...)
    } else {
        antisense <- metaGene(features=fRev, plusCVG=plusCVG, minusCVG=minusCVG,
            size=size, up=up, down=down, ...)
        antisense <- antisense/length(features)
    }
    message("OK")
    
    sense <- sense*normCounts
    antisense <- antisense*normCounts
    return(list(sense=sense, antisense=antisense))
}


samplingMetaGene <- function(features, plusCVG, minusCVG, size=100L, up=10000L,
    down=NULL, nSampling=1000L, samplingRatio=0.1, ...) {
    samplingSize <- round(length(features)*samplingRatio)

    metaList <- mclapply(1:length(features), function(x) {
        metaGene(features=features[x,], plusCVG=plusCVG, minusCVG=minusCVG, 
        size=size, up=up, down=down)
    }, ...)

    allSamples <- mclapply(1:nSampling, function(x) {
        inx <- sample(1:length(features), size=samplingSize, replace=TRUE)
        onesample <- metaList[inx]
        mat <- sapply(onesample, function(x) as.integer(x))
        Rle(apply(mat, 1, sum))
    }, ...)

    M <- sapply(allSamples, function(x) as.integer(x))
    return(Rle(apply(M, 1, median)) / samplingSize)
}


#' Returns a matrix, with rows representing read counts across a specified 
#' gene, or other features of interest.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors, use the argument 'mc.cores'.
#'
#' @param features A GRanges object representing a set of genomic coordinates. 
#' @param reads A GRanges object representing a set of mapped reads. 
#' @param size The size of the moving window.
#' @param up Distance upstream of each f to align and histogram Default: 1 kb.
#' @param down Distance downstream of each f to align and histogram 
#' Default: same as up.
#' @param debug If set to TRUE, provides additional print options. 
#' Default: FALSE
#' @param ... Extra argument passed to mclapply
#' @return Returns a vector representing the 'typical' signal across 
#' genes of different length.
#' @author Charles G. Danko and Minho Chae
##  Returns a matrix of counts.  Rows represent different streches of DNA.
##  Columns represent positions relative to a certain feature.  Summed together,
##  these should be a meta-gene.
##
##  Arguments:
##  f   -> data.frame of: CHR, START, STRAND.
##  p   -> data.frame of: CHR, START, END, STRAND.
##  size    -> The size of the moving window.
##  up  -> Distance upstream of each f to align and histogram.
##  down    -> Distance downstream of each f to align and histogram (NULL).
##
##  Assumptions: Same as MetaGene
metaGeneMatrix <- function(features, reads, size= 50, up=1000, down=up, 
    debug=FALSE, ...) {

    C <- sort(unique(as.character(seqnames(features))))

    ## Run parallel version.
    mcp <- mclapply(seq_along(C), metaGeneMatrix_foreachChrom, C=C, 
        features=features, reads=reads, size=size, up=up, down=down, 
        debug=debug, ...)
    
    ## Append data from all chromosomes.
    H <- NULL
    for(i in seq_along(C)) {
        # Which KG?  prb?
        indxF   <- which(as.character(seqnames(features)) == C[i])
        indxPrb <- which(as.character(seqnames(reads)) == C[i])

        if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
            H <- rbind(H, mcp[[i]])
        }
    }

    return(H)
}


metaGeneMatrix_foreachChrom <- function(i, C, features, reads, size, up, down, 
    debug) {
        # Which KG?  prb?
        indxF   <- which(as.character(seqnames(features)) == C[i])
        indxPrb <- which(as.character(seqnames(reads)) == C[i])

        if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
            # Order -- Make sure, b/c this is one of our main assumptions.  
            # Otherwise violated for DBTSS.
            ord <- order(start(features[indxF,]))
            # Type coersions.
            FeatureStart    <- start(features[indxF,][ord])
            FeatureStr  <- as.character(strand(features[indxF,][ord]))
            PROBEStart  <- start(reads[indxPrb,])
            PROBEEnd    <- end(reads[indxPrb,])
            PROBEStr    <- as.character(strand(reads[indxPrb,]))
            size        <- as.integer(size)
            up          <- as.integer(up)
            down        <- as.integer(down)

            # Set dimensions.
            dim(FeatureStart)   <- c(NROW(FeatureStart), NCOL(FeatureStart))
            dim(FeatureStr)     <- c(NROW(FeatureStr),   NCOL(FeatureStr))
            dim(PROBEStart)     <- c(NROW(PROBEStart),   NCOL(PROBEStart))
            dim(PROBEEnd)       <- c(NROW(PROBEEnd),     NCOL(PROBEEnd))
            dim(PROBEStr)       <- c(NROW(PROBEStr),     NCOL(PROBEStr))

            if(debug) {
                message(C[i],": Counting reads in specified region.")
            }
            Hprime <- .Call("MatrixOfReadsByFeature", FeatureStart, FeatureStr, 
                            PROBEStart, PROBEEnd, PROBEStr, 
                            size, up, down, PACKAGE = "groHMM")
            return(Hprime)
        }
    return(integer(0))
}



#' Returns a histogram of the number of reads in each section of a moving 
#' window of #' variable size across genes.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors, use the argument 'mc.cores'.
#'
#' @param features A GRanges object representing a set of genomic coordinates. 
#' @param reads A GRanges object representing a set of mapped reads. 
#' @param n_windows The number of windows to break genes into.
#' @param debug If set to TRUE, provides additional print options. 
#' Default: FALSE
#' @param ... Extra argument passed to mclapply
#' @return Returns a vector representing the 'typical' signal across genes of 
#' different length.
#' @author Charles G. Danko and Minho Chae
##  Returns a histogram of the number of reads in each section of a
##  moving window of variable size across genes.
##
##  Arguments:
##  f   -> data.frame of: CHR, START, END, STRAND.
##  p   -> data.frame of: CHR, START, END, STRAND.
##  n_windows   -> The resolution of the MetaGene -- i.e. the number of moving 
##  windows to break it into..
##
##  Assumptions:
##  (1) Gene list should be ordered!  
##  (2) Gene list should be pretty short, as most of the processing and 
##  looping over genes is currently done in R.
#
metaGene_nL <- function(features, reads, n_windows=1000, debug=FALSE, ...) {
    C <- sort(unique(as.character(seqnames(features))))
    H <- rep(0,n_windows)
    for(i in 1:NROW(C)) {
        if(debug) {
            message(C[i])
        }

        # Which KG?  prb?
        indxF   <- which(as.character(seqnames(features)) == C[i])
        indxPrb <- which(as.character(seqnames(reads)) == C[i])

        if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
            # Order -- Make sure, b/c this is one of our main assumptions.  
            # Otherwise violated for DBTSS.
            ord <- order(start(features[indxF,])) 

            # Type coersions.
            FeatureStart    <- start(features[indxF,][ord])
            FeatureEnd  <- end(features[indxF,][ord])
            FeatureStr  <- as.character(strand(features[indxF,][ord]))
            PROBEStart  <- start(reads[indxPrb,])
            PROBEEnd    <- end(reads[indxPrb,])
            PROBEStr    <- as.character(strand(reads[indxPrb,]))

            # Set dimensions.
            dim(FeatureStart)   <- c(NROW(FeatureStart), NCOL(FeatureStart))
            dim(FeatureStr)     <- c(NROW(FeatureStr),   NCOL(FeatureStr))
            dim(PROBEStart)     <- c(NROW(PROBEStart),   NCOL(PROBEStart))
            dim(PROBEEnd)       <- c(NROW(PROBEEnd),     NCOL(PROBEEnd))
            dim(PROBEStr)       <- c(NROW(PROBEStr),     NCOL(PROBEStr))

            #for(iFeatures in 1:NROW(FeatureStart)) {
            mcpg <- mclapply(c(1:NROW(FeatureStart)), function(iFeatures) {
                ws <- (FeatureEnd[iFeatures]-FeatureStart[iFeatures])/n_windows 
                ## This WILL be an interger.
                if(debug) {
                    message(C[i],": Counting reads in specified region:",
                            FeatureStart[iFeatures],"-",FeatureEnd[iFeatures])
                    message(C[i],": Window size:",
                            FeatureStart[iFeatures],"-",FeatureEnd[iFeatures])
                    message(C[i],": End-Start:",
                            FeatureEnd[iFeatures]-FeatureStart[iFeatures])
                }
                DataByOne <- .Call("WindowAnalysis", PROBEStart, PROBEEnd, 
                        PROBEStr, FeatureStr[iFeatures], as.integer(1), 
                        as.integer(1), FeatureStart[iFeatures], 
                        FeatureEnd[iFeatures], 
                        PACKAGE = "groHMM")

                if(debug) {
                    message("DataByOne size:",NROW(DataByOne))
                }

                ## This seems almost immeidate on my pentium M machine.
                Hprime <- unlist(lapply(1:NROW(H), function(i) {
                    indx <- ceiling(ws*(i-1)+1):ceiling(ws*i)
                    return(sum(DataByOne[indx]))
                }))

                ## Reverse it for "-" strands.
                if(FeatureStr[iFeatures] == "-")
                    Hprime <- rev(Hprime)

                #H <- H + Hprime/ws ## Put in units of: number of reads/base
                return(Hprime/ws)
            }, ...)

            ## Add genes from mclapply together.
            for(iFeatures in 1:NROW(FeatureStart)) { 
                H<- H+mcpg[[i]] 
            }
        }
        if(debug) {
            message(C[i],": Done!")
        }
    }

    return(H)
}


#' Returns the average profile of tiling array probe intensity values or 
#' wiggle-like count data centered on a set of genomic positions 
#' (specified by 'Peaks').
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors, use the argument 'mc.cores'.
#'
#' @param ProbeData Data.frame representing chromosome, window center, 
#' and a value.
#' @param Peaks Data.frame representing chromosome, and window center.
#' @param size Numeric.  The size of the moving window. Default: 50 bp.
#' @param bins The bins of the meta gene -- i.e. the number of moving windows 
#' to break it into. Default +/- 1kb from center.
#' @return A vector representing the 'typical' signal centered on the peaks of 
#' interest.
#' @author Charles G. Danko and Minho Chae
##  Returns the average profile of tiling array probe intensity values or 
## wiggle-like count data centered on a set of genomic positions.
##
##  Arguments:
##  Peaks       -> data.frame of: CHR, CENTER, STRAND. 
##      (note that STRAND is currenly not supported, and does nothing).
##  ProbeData   -> data.frame of: CHR, CENTER, VALUE
##  bins        -> The bins of the meta gene -- i.e. the number of 
##          moving windows to break it into.
##
##  TODO: 
##  (1) Implement support for a Peaks$starnd 
##  (2) ...
averagePlot <- function(ProbeData, Peaks, size=50, bins= seq(-1000,1000,size)){

    ## For each chromsome.  
    ProbeData$minDist <- rep(999)
    for(chr in unique(Peaks[[1]])) {

      ##   Get the set of data that fits.  
      indxPeaks <- which(Peaks[[1]] == chr)
      
      ## The '$' ensures that chrom is end of line.  Otherwise, grep for chr1 
      ## returns chr10-chr19 as well.
      ## Should work fine, even when chromsome is simply "chrN".
      indxAffxProbes <- grep(paste(chr,"$", sep=""), ProbeData[[1]], perl=TRUE)

      ##   Calculate the minimum distance between the probe and the vector 
      ## over all features ... 
      ProbeData$minDist[indxAffxProbes] <- 
        unlist(lapply(indxAffxProbes, function(x) {
            ## TODO: For strand to switch it, just multiply by strand here.
            return((ProbeData[x,2] - 
                Peaks[indxPeaks,2])[which.min(abs(ProbeData[x,2] - 
                Peaks[indxPeaks,2]))])}))
    }

    ## Make bins.  Take averages and all that...
    means <- unlist(lapply(c(1:NROW(bins)), 
        function(i){mean(ProbeData[(ProbeData$minDist >= 
            (bins[i]-size) & ProbeData$minDist < bins[i]),3])}))
    return(data.frame(windowCenter= bins+(size/2), means))
}
