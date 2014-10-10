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

#' getTxDensity Calculates transcript density.
#'
#' Calculates transcript density for transcripts which overlapps with 
#' annotations.  
#' For 'run genes together' or 'broken up a single annotation' errors, 
#' best overlapped transcripts or annotations are used.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors
#' set the option 'mc.cores'.
#'
#' @param tx GRanges of transcripts. 
#' @param annox GRanges of non-overlapping annotatoins.
#' @param plot Logical.  If TRUE, plot transcript density.  Default: TRUE
#' @param scale Numeric. Scaled size of a gene for transcript density 
#' calculation. 
#' Default: 1000L
#' @param nSampling Numeric. Number of subsampling.  Default: 0L
#' @param samplingRatio Numeric. Ratio of sampling for annotations.  
#' Default: 0.1
#' @param ... Extra argument passed to mclapply.
#' @return Returns a list of FTD, TTD, PostTTS, and AUC. 
#' @author Minho Chae
#' @examples
#' tx <- GRanges("chr7", IRanges(start=seq(1000,4000, by=1000), 
#' width=seq(1000, 1300, by=100)), strand=rep("+", 4))
#' annox <- GRanges("chr7", IRanges(start=seq(1100,4100, by=1000), 
#' width=seq(900, 1200, by=100)), strand=rep("+", 4))
#' ## Not run:
#' # density <- getTxDensity(tx, annox) 
getTxDensity <- function(tx, annox, plot=TRUE, scale=1000L, nSampling=0L, 
    samplingRatio=0.1, ...) {
    ol <- findOverlaps(tx, annox)

    # Count tx
    runGenes <- unique(queryHits(ol[duplicated(queryHits(ol)),]))       
    # Count annox
    brokenUp <- unique(subjectHits(ol[duplicated(subjectHits(ol)),]))   
    cat("Merged annotations: ", length(runGenes), "\n")
    cat("Dissociated a single annotation: ", length(brokenUp), "\n")
    cat("Overlaps between transcript and annotation:", "\n")
    cat("Total = ", length(ol))

    # For each annox, find the best matching tx, runGenes case...
    intx_rg <- pintersect(tx[queryHits(ol),], annox[subjectHits(ol),])
    intx_rg_df <- data.frame(tx=queryHits(ol), annox=subjectHits(ol),   
            oRatio=width(intx_rg)/width(tx[queryHits(ol),]))

    #Tx matches multiple times on a same annox
    remove_rg <- unique(unlist(lapply(runGenes, function(x) {
                    inx <- which(queryHits(ol) == x)
                    m <- which.max(intx_rg_df$oRatio[inx])
                    inx[-m]
                })))

    if (length(remove_rg) > 0)  
        ol <- ol[-remove_rg,]

    # For each transcript, find the best matching annox, brokenUp case..., 
    brokenUp <- unique(subjectHits(ol[duplicated(subjectHits(ol)),]))

    intx_bu  <- pintersect(tx[queryHits(ol),], annox[subjectHits(ol),])
    intx_bu_df <- data.frame(tx=queryHits(ol), annox=subjectHits(ol),   
            oRatio=width(intx_bu)/width(annox[subjectHits(ol),]))


    # Annox matches multiple times on a same transcript
    remove_bu <- unique(unlist(lapply(brokenUp, function(x) {
                    inx <- which(subjectHits(ol) == x)
                    m <- which.max(intx_bu_df$oRatio[inx])
                    inx[-m]
                })))
    if (length(remove_bu) > 0)  
        ol <- ol[-remove_bu,]

    cat(" Used for density = ", length(ol), "\n")

    olTx <- tx[queryHits(ol),]
    # Now get the coverage of selected transcripts
    olStrand <- as.character(strand(tx[queryHits(ol),]))
    olChrom <- seqnames(tx[queryHits(ol),])

    # Get the extended region for annox
    up <- 1L
    down <- 2L
    message("Calculate overlapping ... ", appendLF=FALSE)
    promo  <- unlist(GRangesList(mclapply(subjectHits(ol), function(x) {
                        w <- width(annox[x,])
                        promoters(annox[x,], upstream=round(w*up), 
                            downstream=round(w*down))
                    }, ... )))

    pintx <- pintersect(promo, olTx)

    # Get theoverlapped coverage
    olcvg <- mclapply(1:length(ol), function(x) {
                    t <- olTx[x,]
                    p <- promo[x,]
                    i <- pintx[x,]
                    # Position is relative to the minimum start 
                    minStart <- min(start(t), start(p))
                    t <- shift(t, -minStart+1)
                    p <- shift(p, -minStart+1)
                    i <- shift(i, -minStart+1)
                    r <- reduce(c(t, p, ignore.mcols=TRUE))
                    rTF <- logical(length=width(r))
                    rTF[start(i):end(i)] <- TRUE
                    if (olStrand[x] == "+")
                        Rle(rTF[start(p):end(p)])
                    else
                        rev(Rle(rTF[start(p):end(p)]))
            }, ...)
    message("OK")
                
    message("Scale overlapping ... ", appendLF=FALSE)
    # Get the scaled coverage 
    cvgWidth <- round(up*scale) + round(down*scale)
    sccvg <- mclapply(olcvg, function(x) {
                        getLIValues(x, cvgWidth)
                    }, ... )
    message("OK")

    M <- sapply(sccvg, function(x) as.integer(x))
    sSize <- round(length(ol)*samplingRatio)
    if (nSampling > 0) {
        message("Sampling ... ", appendLF=FALSE)
        allSamples <- mclapply(1:nSampling, function(x) {
            inx <- sample(1:length(ol), size=sSize, replace=TRUE)
            onesample <- M[,inx]
            Rle(apply(onesample, 1, sum))
        }, ...)
        mat <- sapply(allSamples, function(x) as.integer(x))
        profile <- apply(mat, 1, mean)/sSize
        message("OK")
    } else {
        profile <- apply(M, 1, sum)/length(ol)
    }
    if (plot) {
        plot(-(up*scale):(down*scale-1), profile, ylim=c(0, 1), type="l", 
            xlab="Relative to TSS", ylab="Density")
        abline(v=0, col="blue", lty=2)
        abline(v=up*scale, col="blue", lty=2)
    }

    trap.rule <- function(x,f) {sum(diff(x)*(f[-1]+f[-length(f)]))/2}
    # Thanks to: http://tolstoy.newcastle.edu.au/R/help/05/08/9625.html

    FTD <- trap.rule(1:scale, profile[1:scale])/scale
    TTD <- trap.rule(1:scale, profile[(scale+1):(scale*2)])/scale
    PostTTS <- trap.rule(1:scale, profile[(scale*2+1):(scale*3)])/scale
    AUC <- (TTD + (TTD - TTD*FTD))/(1 + TTD)
    #message("FTD: ", round(FTD, 2), " TTD: ", round(TTD, 2), " PostTTS: ", 
    #           round(PostTTS, 2), " AUC: ", round(AUC, 2))
    return(list(FTD=FTD, TTD=TTD, PostTTS=PostTTS, AUC=AUC))
}

getWP <- function (lv, lw) {
    return( ((lv-1)/(lw-1)) * (0:(lw-1)) + 1)
}

getLIValue <- function (x0, y0, x1, y1, x) {
    alpha = (x - x0) / (x1 - x0)
    return(y0 + alpha * (y1 - y0))
}

# Get linear interpolation data
getLIValues <- function (vals, n) {
    vals <- as.integer(vals)
    wp <- getWP(length(vals), n)
    #print(wp)
    result <- seq(n)
    if (length(vals) == 1) {
        result[1:length(result)] <- vals
        return(result)
    }

    result[1] <- vals[1]
    result[length(result)] <- vals[length(vals)]

    if (n > 2) {
        for(i in 2:(n-1)) {
            x <- wp[i]
            x0 <- floor(x)
            result[i] <- getLIValue(x0, vals[x0], x0+1, vals[x0+1], x)
        }
    }
    return(Rle(round(result)))
}



#' evaluateHMM Evaluates HMM calling. 
#'
#' Evaluates HMM calling of transripts compared to known annotations. 
#'
#' @param tx GRanges of transcripts predicted by HMM. 
#' @param annox GRanges of non-overlapping annotatoins.
#' @return a list of error information; merged annotations, dissociated annotation, 
#' total, and rate.
#' @author Minho Chae
#' @examples
#' tx <- GRanges("chr7", IRanges(start=seq(100, 1000, by=200), 
#' width=seq(100, 1000, by=100)), strand="+")
#' annox <- GRanges("chr7", IRanges(start=seq(110, 1100, by=150), 
#' width=seq(100, 1000, by=150)), strand="+")
#' error <- evaluateHMMInAnnotations(tx, annox)
evaluateHMMInAnnotations <- function (tx, annox) {
    o <- findOverlaps(tx, annox)
    # count tx
    merged <- length(unique(queryHits(o[duplicated(queryHits(o)),])))  
    # count annox
    dissociated <- length(unique(subjectHits(o[duplicated(subjectHits(o)),]))) 

    eval <- data.frame(merged=merged, dissociated=dissociated,
            total=(merged+ dissociated), errorRate=(merged + dissociated)/
                (length(tx) + length(annox)),
            txSize =length(tx))

    intx <- pintersect(tx[queryHits(o),], annox[subjectHits(o),])
    overlap <- data.frame(txRatio=width(intx)/width(tx[queryHits(o),]),
                annoxRatio= width(intx)/width(annox[subjectHits(o),]),
                overBases=width(intx), 
                similarity=width(intx)/pmax(width(tx[queryHits(o),]), 
               width(annox[subjectHits(o),])))

    return(list(eval=eval, overlap=overlap))
}
