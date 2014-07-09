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
plot2Ranges <- function(tr, gr, main = "NA", col = "black", sep = 0.5, ...) {
    height <- 1
    xlim <- c(min(min(start(tr), min(start(gr)))), max(max(end(tr)), 
                max(end(gr)))) 

    plot.new()
    ybottom <- seq(1,length(tr)+length(gr))
    plot.window(xlim, c(0,max(ybottom+height)))

    if (length(tr)>0) {
        tbottom <- seq(1,length(tr))
        rect(start(tr)-sep, tbottom, end(tr)+sep, tbottom + height-sep, 
                col = "gray", ...)
    }

    bSymbol <- FALSE
    if ("symbol" %in% colnames(elementMetadata(gr)))
    bSymbol <- TRUE

    if (length(gr) >0) {
        gbottom <- seq(1,length(gr))
        rect(start(gr)-sep, gbottom, end(gr)+sep, gbottom + height-sep,
            col = rgb(red=255, green=0, blue=0, alpha=150, maxColorValue=255))
        if (bSymbol)
            text((start(gr) + end(gr))/2, gbottom + 0.2, 
                elementMetadata(gr)$symbol) 
    }
    title(main)
    axis(1)
}
    


#
# Break an Interval by break points with minimum gap 
#------------------------------------------------------------------------------
breakInterval <- function(gr, brPos, gap=5, strand="+") {
    result <- rep(gr, length(brPos)+1)
    if (strand == "+") {
        for(i in seq_along(brPos)) { 
            end(result[i,]) <- brPos[i] - gap   
            start(result[i+1,]) <- brPos[i]
        }
    } else {
        for(i in seq_along(brPos)) { 
            end(result[i,]) <- brPos[i] 
            start(result[i+1,]) <- brPos[i] + gap
        }
    }
    return(result)
}


#' breakTranscriptsOnGenes Breaks transcripts on genes 
#'
#' Breaks transcripts when they are overlapped with multiple well annotated 
#' genes.
#'
#' @param tx GRanges of transcripts.
#' @param annox GRanges of non-overlapping annotations for reference.
#' @param strand Takes "+" or "-" Default: "+"
#' @param geneSize Numeric. Minimum gene size in annox to be used as reference. 
#' Default: 5000
#' @param threshold Numeric. Ratio of overlapped region relative to a gene 
#' width. 
#' Transcripts only greater than this threshold are subjected to be broken. 
#' Default: 0.8
#' @param gap Numeric.  Gap (bp) between broken transcripts.  Default: 5
#' @param plot Logical.  If set to TRUE, show each step in a plot. 
#' Default: FALSE
#' @author Minho Chae and Charles G. Danko
#' @return Returns GRanges object of broken transcripts. 
#' @examples
#' tx <- GRanges("chr7", IRanges(1000, 30000), strand="+")
#' annox <- GRanges("chr7", IRanges(start=c(1000, 20000), 
#'              width=c(10000,10000)), strand="+")
#' bPlus <- breakTranscriptsOnGenes(tx, annox, strand="+")
breakTranscriptsOnGenes <- function(tx, annox, strand="+", geneSize=5000, 
    threshold=0.8, gap=5, plot=FALSE) {
    tx <- tx[strand(tx) == strand,]
    annox <- annox[strand(annox) == strand,]
    mcols(tx)$status <- "NA"

    annox <- annox[width(annox) > geneSize,]
    ol <- findOverlaps(tx, annox)

    ol.df <- data.frame(trans=queryHits(ol), gene=subjectHits(ol), 
        ratio=width(pintersect(tx[queryHits(ol),], annox[subjectHits(ol),]))
        /width(annox[subjectHits(ol),]))

    # Keep only major overlappings
    ol.df <- ol.df[ol.df$ratio > threshold,]

    # Keep only multiple genes on a transcript cases
    dupTrans <- unique(ol.df$trans[duplicated(ol.df$trans)])
    ol.df <- ol.df[ol.df$trans %in% dupTrans,]

    ol.tab <- table(ol.df$trans)
    ol.tabCS <- cumsum(ol.tab)

    bT <- NULL

    for (i in seq_along(ol.tab)) {
        txNo <- as.integer(names(ol.tab[i]))
        aNo <- ol.df$gene[(ol.tabCS[i]-ol.tab[i]+1):ol.tabCS[i]]

        if (strand == "+") {
            # No need of break position for the first annotation
            brPos <- start(annox[aNo[-1],])        
        } else {
            # No need of break position for the last annotation
            brPos <- end(annox[aNo[-length(aNo)],])  
        }

        frags <- breakInterval(tx[txNo,], brPos=brPos, gap=gap, strand=strand)

        if(is.null(bT)) {
            bT <- frags 
        } else {
            bT <- c(bT, frags) 
        }
        if (plot) {
            par(mfrow=c(2,1))
            plot2Ranges(tx[txNo,], annox[aNo,], main="Before")
            plot2Ranges(frags, annox[aNo,], main="After")
            Sys.sleep(5)
        }
    }

    mcols(bT)$status <- "broken"
    mcols(bT)$ID <- paste(seqnames(bT), "_", start(bT), strand(bT), sep="")
    okTrans <- setdiff(1:length(tx), dupTrans)
    all <- c(tx[okTrans,], bT)
    cat(length(unique(ol.df$trans)), " transcripts are broken into ", 
        length(bT), "\n")

    return(all[order(as.character(seqnames(all)), start(all)),])
}

#' combineTranscripts Combines transnscipts. 
#'
#' Combines transcripts  that are within the same gene annotation, combining 
#' smaller transcripts for genes
#'  with low regulation into a single transcript representing the gene.
#'
#' @param tx GRanges of transcripts.
#' @param annox GRanges of non-overlapping annotations for reference.
#' @param geneSize Numeric. Minimum gene size in annotations to be used as 
#' reference. 
#' Default: 1000
#' @param threshold Numeric. Ratio of overlapped region relative to transcript
#' width. 
#' Transcripts only greater than this threshold are subjected to be combined. 
#' Default: 0.8
#' @param plot Logical.  If set to TRUE, show easch step in a plot. 
#' Default: FALSE
#' @return Returns GRanges object of combined transcripts. 
#' @author Minho Chae and Charles G. Danko
#' @examples
#' tx <- GRanges("chr7", IRanges(start=c(1000, 20000), width=c(10000,10000)), 
#' strand="+")
#' annox <- GRanges("chr7", IRanges(1000, 30000), strand="+")
#' combined <- combineTranscripts(tx, annox)
combineTranscripts <- function(tx, annox, geneSize=1000, threshold=0.8, 
    plot=FALSE) {
    annox <- annox[width(annox) > geneSize,]
    ol <- findOverlaps(tx, annox)
    ol.df <- data.frame(trans=queryHits(ol), gene=subjectHits(ol),
        ratio=width(pintersect(tx[queryHits(ol),], annox[subjectHits(ol),]))
                /width(tx[queryHits(ol),]))

    # Keep only major overlappings
    ol.df <- ol.df[ol.df$ratio > threshold,]

    # Keep only multiple transcripts on a gene cases
    dupGenes <- ol.df$gene[which(duplicated(ol.df$gene))]
    ol.df <- ol.df[ol.df$gene %in% dupGenes,]

    uniqGene <- unique(ol.df$gene)
    N <- length(uniqGene)
    cT <- GRanges(seqnames=Rle(rep("chr1", N)),
        ranges = IRanges(rep(1, N), rep(max(end(tx)), N)),
        strand = rep("+", N),
        type = Rle(rep("tx", N)))
    seqlevels(cT) <- seqlevels(tx)

    for (i in seq_along(uniqGene)) {
        block <- ol.df[uniqGene[i]==ol.df$gene,]
        strand(cT[i,]) <- strand(tx[block[1,"trans"],])
        seqnames(cT[i,]) <- seqnames(tx[block[1,"trans"],])
        start(cT[i,]) <- start(tx[block[1,"trans"],])
        end(cT[i,]) <- end(tx[block[NROW(block),"trans"],])
        if (plot) { # plot before and after the comine
            #b4combine <- c(gr[uniqGene[i],], tr[ol.df[ol.df$gene == 
            # uniqGene[i], "trans"],-2])
            #af.combine <- c(gr[uniqGene[i],], cT[i,])
            par(mfrow=c(2,1))
            plot2Ranges(tx[ol.df[ol.df$gene == uniqGene[i], "trans"],], 
                annox[uniqGene[i],], main="Before")
            plot2Ranges(cT[i,], annox[uniqGene[i],], main="After")
            Sys.sleep(5)
        }
    }

    cat(NROW(ol.df), " transcripts are combined to ", NROW(cT), "\n")
    mcols(cT)$ID <- paste(seqnames(cT), "_", start(cT), strand(cT), sep="")
    mcols(cT)$status <- "combined"

    okTrans <- setdiff(seq_along(tx), ol.df$trans)
    all <- c(tx[okTrans,], cT)
    return(all[order(as.character(seqnames(all)), start(all)),])
}
