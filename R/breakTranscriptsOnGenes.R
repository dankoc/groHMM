###########################################################################
##
##   Copyright 2012, 2013 Minho Chae.
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
	xlim <- c(min(min(start(tr), min(start(gr)))), max(max(end(tr)), max(end(gr))))

	plot.new()
	ybottom <- seq(1,length(tr)+length(gr))
	plot.window(xlim, c(0,max(ybottom+height)))

	if (length(tr)>0) {
		tbottom <- seq(1,length(tr))
		rect(start(tr)-sep, tbottom, end(tr)+sep, tbottom + height-sep, col = "gray", ...)
	}

	bSymbol <- FALSE
	if ("symbol" %in% colnames(elementMetadata(gr)))
	bSymbol <- TRUE

	if (length(gr) >0) {
		gbottom <- seq(1,length(gr))
		rect(start(gr)-sep, gbottom, end(gr)+sep, gbottom + height-sep,
			col = rgb(red=255, green=0, blue=0, alpha=150, maxColorValue=255))
		if (bSymbol)
			text((start(gr) + end(gr))/2, gbottom + 0.2, elementMetadata(gr)$symbol)
	}
	title(main)
	axis(1)
}
	


#
# Break an Interval by break points with minimum gap 
#----------------------------------------------------------------------------------------------
breakInterval <- function(gr, brPos, gap=5, strand="+") {
	result <- rep(gr, NROW(brPos)+1)
	if (strand == "+") {
		for(i in 1:NROW(brPos)) { 
		    end(result[i,]) <- brPos[i] - gap   
		    start(result[i+1,]) <- brPos[i]
		}
	} else {
		for(i in 1:NROW(brPos)) { 
		    end(result[i,]) <- brPos[i] 
		    start(result[i+1,]) <- brPos[i] + gap
		}
	}
	return(result)
}

# Break Transcrips on Genes
#----------------------------------------------------------------------------------------------
#' breakTranscriptsOnGenes Breaks transcripts on genes 
#'
#' Breaks transcripts when they are fully overlap multiple well annotated genes.
#'
#' @param tx GRanges of transcripts.
#' @param annotations GRanges of non-overlapping annotations for reference.
#' @param strand Takes "+" or "-" Default: "+"
#' @param geneSize Numeric. Minimum gene size in gr to be used as reference. Default: 5000
#' @param threshold Numeric. Threshold for calling the gene part of the transcript.  Default: 0.8
#' @param gap Numeric.  Gap (bp) between broken transcripts.  Default: 5
#' @param debug Logical.  If set to TRUE, show easch step in a plot. Default: FALSE
#' @author Minho Chae and Charles G. Danko
breakTranscriptsOnGenes <- function(tx, annotations, strand="+", geneSize=5000, threshold=0.8, gap=5, debug=FALSE) {
	tr <- tx
	gr <- annotations
	tr <- tr[as.character(strand(tr)) == strand,]
	gr <- gr[as.character(strand(gr)) == strand,]
	elementMetadata(tr)$fixError <- "NA"

	print(paste("Initial transcripts:", length(tr)))

	gr <- gr[width(gr) > geneSize,]
	ol <- findOverlaps(tr, gr)

	ol.df <- data.frame(trans=queryHits(ol), gene=subjectHits(ol), 
	    geneWidthRate = width(pintersect(tr[queryHits(ol),], gr[subjectHits(ol),]))
		/width(gr[subjectHits(ol),]))

	# remove by threshold
	ol.df <- ol.df[ol.df$geneWidthRate > threshold,]

	# keep only multiple genes on a transcript cases
	dupTrans <- unique(ol.df$trans[duplicated(ol.df$trans)])
	ol.df <- ol.df[ol.df$trans %in% dupTrans,]

	ol.tab <- table(ol.df$trans)
	ol.tabCS <- cumsum(ol.tab)

	bT <- NULL

	for (i in 1:NROW(ol.tab)) {
		txNo <- as.integer(names(ol.tab[i]))
		gNo <- ol.df$gene[(ol.tabCS[i]-ol.tab[i]+1):ol.tabCS[i]]

		if (strand == "+") {
			brPos <- start(gr[gNo[-1],])        # No need of break position for the first annotation
		} else {
			brPos <- end(gr[gNo[-NROW(gNo)],])  # No need of break position for the last annotation
		}

		frags <- breakInterval(tr[txNo,], brPos=brPos, gap=gap, strand=strand)

		if(is.null(bT)) {
			bT <- frags 
		} else {
			bT <- c(bT, frags) 
		}
		if (debug) {
			par(mfrow=c(2,1))
			plot2Ranges(tr[txNo,], gr[gNo,], main="Before")
			plot2Ranges(frags, gr[gNo,], main="After")
			browser()
		}
	}

	elementMetadata(bT)$fixError <- "broken"
	elementMetadata(bT)$ID <- paste(seqnames(bT), "_", start(bT), strand(bT), sep="")
	okTrans <- setdiff(1:NROW(tr), dupTrans)
	all <- c(tr[okTrans,], bT)
	print(paste(NROW(unique(ol.df$trans)), "transcripts are broken into", length(bT)))
	print(paste("Final transcripts:", length(all)))

	return(all[order(as.character(seqnames(all)), start(all)),])
}

#-----------------------------------------------------------------------------------------------
# Combine Transcrips on Genes
# To do: fix the name of the combined transcripts
#-----------------------------------------------------------------------------------------------
#' combineTranscripts Combines transnscipts. 
#'
#' Combines transcripts  that are within the same gene annotation, combining smaller transcripts for genes
#'  with low regulation into a single transcript representing the gene.
#'
#' @param tx GRanges of transcripts.
#' @param annotations GRanges of non-overlapping annotations for reference.
#' @param geneSize Numeric. Minimum gene size in annotations to be used as reference. Default: 1000
#' @param threshold Numeric. Threshold for calling the gene part of the transcript.  Default: 0.8
#' @param debug Logical.  If set to TRUE, show easch step in a plot. Default: FALSE
#' @author Minho Chae and Charles G. Danko
combineTranscripts <- function(tx, annotations, geneSize=1000, threshold=0.8, debug=FALSE) {
	tr <- tx
	gr <- annotations
	print(paste("Initial transcripts:", length(tr)))
	gr <- gr[width(gr) > geneSize,]
	ol <- findOverlaps(tr, gr)
	ol.df <- data.frame(trans=queryHits(ol), gene=subjectHits(ol),
	    trWidthRate = width(pintersect(tr[queryHits(ol),], gr[subjectHits(ol),]))
			    /width(tr[queryHits(ol),]))
	# remove by threshold
	ol.df <- ol.df[ol.df$trWidthRate > threshold,]

	# keep only multiple transcripts on a gene cases
	dupGenes <- ol.df$gene[which(duplicated(ol.df$gene))]
	ol.df <- ol.df[ol.df$gene %in% dupGenes,]

	uniqGene <- unique(ol.df$gene)
	#cT <- gr[uniqGene,-(1:2)]  # template
	#cT <- gr[uniqGene,-1]  # template
	#cT <- # template
	N <- NROW(uniqGene)
	cT <- GRanges(seqnames=Rle(rep("chr1", N)),
	    ranges = IRanges(rep(1, N), rep(max(end(tr)), N)),
	    strand = rep("+", N),
	    type = Rle(rep("tx", N)))
	seqlevels(cT) <- seqlevels(tr)

	for (i in 1:NROW(uniqGene)) {
		block <- ol.df[uniqGene[i]==ol.df$gene,]
		strand(cT[i,]) <- strand(tr[block[1,"trans"],])
		seqnames(cT[i,]) <- seqnames(tr[block[1,"trans"],])
		start(cT[i,]) <- start(tr[block[1,"trans"],])
		end(cT[i,]) <- end(tr[block[NROW(block),"trans"],])
		if (debug) { # plot before and after the comine
			#b4combine <- c(gr[uniqGene[i],], tr[ol.df[ol.df$gene == uniqGene[i], "trans"],-2])
			#af.combine <- c(gr[uniqGene[i],], cT[i,])
			par(mfrow=c(2,1))
			plot2Ranges(tr[ol.df[ol.df$gene == uniqGene[i], "trans"],], gr[uniqGene[i],], main="Before")
			plot2Ranges(cT[i,], gr[uniqGene[i],], main="After")
			browser()
		}
	}

	print(paste(NROW(ol.df), "transcripts are combined to", NROW(cT)))
	elementMetadata(cT)$ID <- paste(seqnames(cT), "_", start(cT), strand(cT), sep="")
	elementMetadata(cT)$fixError <- "combined"

	okTrans <- setdiff(1:NROW(tr), ol.df$trans)
	all <- c(tr[okTrans,], cT)
	print(paste("Final transcripts:", length(all)))
	return(all[order(as.character(seqnames(all)), start(all)),])
}

