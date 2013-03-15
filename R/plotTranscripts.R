
#-----------------------------------------------------------------------------------------------
# plotTranscripts
#-----------------------------------------------------------------------------------------------
#' plotTranscripts Plots transcripts. 
#'
#' Plots transcripts with their associated genes. 
#'
#' @param tx GRanges of transcripts. 
#' @param features GRanges of non-overlapping annotatoins.
#' @param chrom Character.  Target transcript chromosome; NA for all chromosomes.  Default: NA
#' @param strand Character.  "+" or "-";  NA for all chromosomes.  Default: NA
#' @param scale Logical.  If TRUE, plot is scaled to relative to gene size of 30000 bp. Default: TRUE
#' @param runGenes Character.  One of "best", "all", or "none" for the transcripts running over multiple gene annotations. Default: "best"

#' @param brokenAnnotation Character.  One of "best", "all", or "none" for the transcripts breaking one annotation. Defualt: "best"
#' @param first Integer. First n transcripts to plot after transcripts are ordered by their gene sizes. Default: NA
#' @param last Integer. Last n transcripts to plot after transcripts are ordered by their gene sizes. Default: NA
#' @param filename Character.  File name for output. Default: NA
#' @author Minho Chae
plotTranscripts <- function(tx, features, chrom=NA, strand=NA, scale=TRUE, runGenes="best", brokenAnnotation="best", first=NA, last=NA, filename=NA) {
	if (!is.na(chrom)) {
		tx <- tx[as.character(seqnames(tx)) == chrom,]
		features <- features[as.character(seqnames(features)) == chrom,]
	}

	if (is.na(strand)) {  # both strands
		txPlus <- tx[strand(tx)=="+",]
		fePlus <- features[strand(features)=="+",]
		txMinus <- tx[strand(tx)=="-",]
		feMinus <- features[strand(features)=="-",]
		plusResult <- plotTranscriptsStrand(txPlus, fePlus, strand="+", scale=scale, runGenes=runGenes, brokenAnnotation=brokenAnnotation)
		minusResult <- plotTranscriptsStrand(txMinus, feMinus, strand="-", scale=scale, runGenes=runGenes, brokenAnnotation=brokenAnnotation)

		t <- c(plusResult[[1]], minusResult[[1]])
		g <- c(plusResult[[2]], minusResult[[2]])
	} else {
		result <- plotTranscriptsStrand(tx[strand(tx)==strand,], features[strand(features)==strand,],
			strand=strand, scale=scale, runGenes=runGenes, brokenAnnotation=brokenAnnotation)
		t <- result[[1]]
		g <- result[[2]]
	}


	if (!is.na(first)) {
		t <- t[1:first,]
		g <- g[1:first,]
	}

	if (!is.na(last)) {
		t <- t[(NROW(t)-last+1):NROW(t),]
		g <- g[(NROW(g)-last+1):NROW(g),]
	}

	o <- order(width(g))
	t <- t[o,]
	g <- g[o,]
	#browser()
	xlim <- t
	height <- 1
	sep <- 0.5
	if (is(xlim, "GenomicRanges"))
		xlim <- c(min(start(xlim)), max(end(xlim)))

	if (!is.na(filename))
		png(filename)

	plot.new()
	ybottom <- seq(1,length(t))
	plot.window(xlim, c(0, max(ybottom+height)))
	rect(start(t)-sep, ybottom, end(t)+sep, ybottom + height-sep, col = "gray")

	rect(start(g)-sep, ybottom, end(g)+sep, ybottom + height-sep,
	col = rgb(red=255, green=0, blue=0, alpha=150, maxColorValue=255))
	#title(main)
	axis(1)

	if (!is.na(filename))
		dev.off()

	print(paste("No of transcripts:", length(t)))
}

plotTranscriptsStrand <- function(tr, gr, strand, scale, runGenes="best", brokenAnnotation="best") {
	if (length(tr)*length(gr)==0)
		return(NULL)

	ol <- findOverlaps(tr, gr)

	# Run over genes
	dupTrans <- queryHits(ol)[which(duplicated(queryHits(ol)))]
	if (NROW(dupTrans) > 0) {
		print(paste("Run over genes:", NROW(dupTrans)))
		int <- pintersect(tr[queryHits(ol),], gr[subjectHits(ol),])
		#int.df <- cbind(as.data.frame(ol), width(int) / width(tr[queryHits(ol),]))
		int.df <- data.frame(tr=queryHits(ol), gr=subjectHits(ol), oRatio=width(int) / width(tr[queryHits(ol),]))
		#colnames(int.df) <- c("tr", "gr", "oRatio")
		maxOverlaps <- sapply(dupTrans, function (x) {
				    matchInx <- which(int.df$tr == x)
				    maxInx <- which.max(int.df$oRatio[matchInx])
				    matchInx[maxInx] })
		allOverlaps <- which(queryHits(ol) %in% dupTrans)
		if (runGenes == "best") {
		    ol <- ol[-setdiff(allOverlaps, maxOverlaps),]
		} else if (runGenes == "none") {
		    ol <- ol[-allOverlaps,]
		}
	}

	# transcript which broken annotation
	dupGenes <- subjectHits(ol)[which(duplicated(subjectHits(ol)))]
	if (NROW(dupGenes) > 0) {
		print(paste("Broken annotation:", NROW(dupGenes)))
		int <- pintersect(tr[queryHits(ol),], gr[subjectHits(ol),])
		# int.df <- cbind(as.data.frame(ol), width(int) / width(gr[subjectHits(ol),]))
		int.df <- data.frame(tr=queryHits(ol), gr=subjectHits(ol), oRatio=width(int) / width(tr[queryHits(ol),]))
		# colnames(int.df) <- c("tr", "gr", "oRatio")
		maxOverlaps <- sapply(dupGenes, function (x) {
				     matchInx<- which(int.df$gr == x)
				     maxInx <- which.max(int.df$oRatio[matchInx])
				     matchInx[maxInx] })
		allOverlaps <- which(subjectHits(ol) %in% dupGenes)
		if (brokenAnnotation == "best") {
		    ol <- ol[-setdiff(allOverlaps, maxOverlaps),]
		} else if (runGenes == "none") {
		    ol <- ol[-allOverlaps,]
		}
	}

	# t <- tr[as.matrix(ol)[,1],]
	# g <- gr[as.matrix(ol)[,2],]
	t <- tr[queryHits(ol)]
	g <- gr[subjectHits(ol)]

	#write.table(as.data.frame(t), file="transcripts.tmp", quote=FALSE, sep="\t", row.names=TRUE)
	#write.table(as.data.frame(g), file="gene.tmp", quote=FALSE, sep="\t", row.names=TRUE)
	#browser()

	seqlengths(t) <- NA
	seqlengths(g) <- NA
	# make it all relative to gene start
	if (strand == "+") {
		t <- shift(t, -start(g))
		g <- shift(g, -start(g))
	} else {
		t <- shift(t, -end(g))
		g <- shift(g, -end(g))
	}

	if (scale) {
	# let's shift transcript (x) to the scaled position and then resize
		if (strand == "+") {
			trDist <- start(t) - start(t) * (30000/width(g))
			t <- shift(t, -trDist)
		} else {
			trDist <- end(t) - end(t) * (30000/width(g))
			t <- shift(t, -trDist)
		}
		t <- resize(t, width(t)*(30000/width(g)))  # resize anchors at start for +
		g <- resize(g, 30000)  # gene is okay
	}

	if (strand == "-") {  # flip it
		g <- shift(g, -start(g))
		t <- shift(t, -1*(start(t) + end(t)))
	}

	#print(paste("No of transcripts:", length(t)))
	return(list(t, g))
}

#' evaluateHMM Evaluates HMM calling. 
#'
#' Evaluates HMM calling of transripts compared to know annotations. 
#'
#' @param tx GRanges of transcripts predicted by HMM. 
#' @param annotations GRanges of non-overlapping annotatoins.
#' @return  List of evaluation information; runGeneError, brokenError, overlapQuality.
#' @author Minho Chae
evaluateHMM <- function(tx, annotations) {
	gr <- annotations
	ol <- findOverlaps(tx, gr)
	rgError <-  NROW(unique(queryHits(ol[duplicated(queryHits(ol)),])))
	brokenError <- NROW(unique(subjectHits(ol[duplicated(subjectHits(ol)),])))
	print("=======================================")
	print("HMM Errors")
	print("=======================================")
	print(paste("Run genes together:", rgError))
	print(paste("Broken up a single annotation:", brokenError))

	#intx <- pintersect(tx[queryHits(ol),], gr[subjectHits(ol),])
	#oQuality <- data.frame(tID=elementMetadata(tx[queryHits(ol),])$ID,
	#                             gID=elementMetadata(gr[subjectHits(ol),])$ID,
	#                             gSymbol=elementMetadata(gr[subjectHits(ol),])$symbol,
	#                             tOverlap=width(intx)/width(tx[queryHits(ol),]),
	#                             gOverlap=width(intx)/width(gr[subjectHits(ol),]),
	#                             overBases=width(intx),
	#                             similarity=width(intx)/pmax(width(tx[queryHits(ol),]), width(gr[subjectHits(ol),])))

	intx <- pintersect(tx[queryHits(ol),], gr[subjectHits(ol),])
	oQuality <- data.frame(tOverlap=width(intx)/width(tx[queryHits(ol),]),
			    gOverlap=width(intx)/width(gr[subjectHits(ol),]),
			    overBases=width(intx),
			    similarity=width(intx)/pmax(width(tx[queryHits(ol),]), width(gr[subjectHits(ol),])))


	txWidth <- summary(width(tx[queryHits(ol),]))
	grWidth <- summary(width(gr[subjectHits(ol),]))
	txOverlap <- summary(oQuality$tOverlap)
	grOverlap <- summary(oQuality$gOverlap)
	overBases <- summary(oQuality$overBases)
	similarity <- summary(oQuality$similarity)
	print("=======================================")
	print(paste("Overlap quality (Median Mean)"))
	print("=======================================")
	print(paste("Tx length:", txWidth["Median"], txWidth["Mean"]))
	print(paste("Gene length:", grWidth["Median"], grWidth["Mean"]))
	print(paste("Tx overlap:", txOverlap["Median"], txOverlap["Mean"]))
	print(paste("Gene overlap:", grOverlap["Median"], grOverlap["Mean"]))
	print(paste("Overbases:", overBases["Median"], overBases["Mean"]))
	print(paste("Similarity:", similarity["Median"], similarity["Mean"]))

	return(list(runGeneError=rgError, brokenError=brokenError, overlapQuality=oQuality))
}

#-----------------------------------------------------------------------------------------------
# plotTHistogram
#-----------------------------------------------------------------------------------------------
plotTHistogramStrand <- function(tr, gr, strand, scale, runGenes="best", brokenAnnotation="best") {
	# Assume tr and gr has only one strand
	ol <- findOverlaps(tr, gr)
	ini.ol <- ol

	if(length(ol) == 0){
		if (length(tr) == 0) {
			t <- tr
		} else {
			t=tr[-(1:length(tr)),]
		}
		if (length(gr) == 0) {
			g <- gr
		} else {
			g=gr[-(1:length(gr)),]
		}
		return(list(t,g))
	}

	# transcripts which run over genes
	dupTrans <- queryHits(ol)[which(duplicated(queryHits(ol)))]
	if (NROW(dupTrans)>0) {
		print(paste("Runover genes:", NROW(dupTrans)))
		int <- pintersect(tr[queryHits(ol),], gr[subjectHits(ol),])
		#int.df <- cbind(as.data.frame(ol), width(int) / width(tr[queryHits(ol),]))
		int.df <- data.frame(tr=queryHits(ol), gr=subjectHits(ol), oRatio=width(int) / width(tr[queryHits(ol),]))
		#colnames(int.df) <- c("tr", "gr", "oRatio")
		maxOverlaps <- sapply(dupTrans, function (x) {
				     matchInx <- which(int.df$tr == x)
				     maxInx <- which.max(int.df$oRatio[matchInx])
				     matchInx[maxInx] })
		allOverlaps <- which(queryHits(ol) %in% dupTrans)
		if (runGenes == "best") {
			ol <- ol[-setdiff(allOverlaps, maxOverlaps),]
		} else if (runGenes == "none") {
			ol <- ol[-allOverlaps,]
		}
	}

    # transcript which broken annotation
	dupGenes <- subjectHits(ol)[which(duplicated(subjectHits(ol)))]
	if (NROW(dupGenes) > 0) {
		print(paste("Broken annotation:", NROW(dupGenes)))
		int <- pintersect(tr[queryHits(ol),], gr[subjectHits(ol),])
		# int.df <- cbind(as.data.frame(ol), width(int) / width(gr[subjectHits(ol),]))
		int.df <- data.frame(tr=queryHits(ol), gr=subjectHits(ol), oRatio=width(int) / width(tr[queryHits(ol),]))
		# colnames(int.df) <- c("tr", "gr", "oRatio")
		maxOverlaps <- sapply(dupGenes, function (x) {
				     matchInx<- which(int.df$gr == x)
				     maxInx <- which.max(int.df$oRatio[matchInx])
				     matchInx[maxInx] })
		allOverlaps <- which(subjectHits(ol) %in% dupGenes)
		if (brokenAnnotation == "best") {
			ol <- ol[-setdiff(allOverlaps, maxOverlaps),]
		} else if (runGenes == "none") {
			ol <- ol[-allOverlaps,]
		}
	}

	#t <- tr[as.matrix(ol)[,1],]
	#g  <- gr[as.matrix(ol)[,2],]
	t <- tr[queryHits(ol),]
	g  <- gr[subjectHits(ol),]

	o <- order(width(g))
	t <- t[o]
	g <- g[o]

	seqlengths(t) <- NA
	seqlengths(g) <- NA
	# make it all relative to gene start
	if (strand == "+") {
		t <- shift(t, -1*start(g))
		g <- shift(g, -1*start(g))
	} else {
		t <- shift(t, -end(g))
		g <- shift(g, -end(g))
	}

	if (scale) {
		if (strand == "+") {
			trDist <- start(t) - start(t) * (30000/width(g))
			t <- shift(t, -trDist)
		} else {
			trDist <- end(t) - end(t) * (30000/width(g))
			t <- shift(t, -trDist)
		}
		t <- resize(t, width(t)*(30000/width(g)))  # resize anchors at start for +
		g <- resize(g, 30000)  # gene is okay
	}

	if (strand == "-") {  # flip it
		g <- shift(g, -start(g))
		t <- shift(t, -1*(start(t) + end(t)))
	}

	return(list(t, g))
}



#-----------------------------------------------------------------------------------------------
# plotTHistogram
#-----------------------------------------------------------------------------------------------
# tr: txn, gr: annotation
# overlapGenes: c("best", "all", "none")
# overlapTrans: c("best", "all", "none")
#-----------------------------------------------------------------------------------------------
#' plotTranscripts Plots a transcript histogram. 
#'
#'Plots transcripts stacked together for their associated genes.  Gene size can be absolute or scaled to 30000 bp.
#'
#' @param tx GRanges of transcripts. 
#' @param features GRanges of non-overlapping annotatoins.
#' @param chrom Character.  Target transcript chromosome; NA for all chromosomes. Default: NA
#' @param strand Character.  "+" or "-";  NA for all chromosomes.  Default: NA
#' @param scale Logical.  If TRUE, plot is scaled to relative to gene size of 30000 bp. Default: TRUE
#' @param runGenes Character.  One of "best", "all", or "none" for the transcripts running over multiple gene annotations. Default: "best"
#' @param brokenAnnotation Character.  One of "best", "all", or "none" for the transcripts breaking one annotation. Default: "best"
#' @param filename Character.  File name for output. Default: NA
#' @author Minho Chae
plotTHistogram <- function(tx, features, chrom=NA, strand=NA, scale=TRUE, runGenes="best", brokenAnnotation="best", filename=NA) {
	if (!is.na(chrom)) {
		tx <- tx[as.character(seqnames(tx)) == chrom,]
		features <- features[as.character(seqnames(features)) == chrom,]
	}

	if (is.na(strand)) {  # both strands
		txPlus <- tx[strand(tx)=="+",]
		fePlus <- features[strand(features)=="+",]
		txMinus <- tx[strand(tx)=="-",]
		feMinus <- features[strand(features)=="-",]
		plusResult <- plotTHistogramStrand(txPlus, fePlus, strand="+", scale=scale, runGenes=runGenes, brokenAnnotation=brokenAnnotation)
		minusResult <- plotTHistogramStrand(txMinus, feMinus, strand="-", scale=scale, runGenes=runGenes, brokenAnnotation=brokenAnnotation)

		t <- c(plusResult[[1]], minusResult[[1]])
		g <- c(plusResult[[2]], minusResult[[2]])
	} else {
		result <- plotTHistogramStrand(tx[strand(tx)==strand,], features[strand(features)==strand,],
			strand=strand, scale=scale, runGenes=runGenes, brokenAnnotation=brokenAnnotation)
		t <- result[[1]]
		g <- result[[2]]
	}

	trap.rule <- function(x,f) {sum(diff(x)*(f[-1]+f[-length(f)]))/2}
	# source: http://tolstoy.newcastle.edu.au/R/help/05/08/9625.html

	hist <- numeric(max(end(t)- min(start(t)))+1)
	#histMat <- matrix(0, nrow=length(g), ncol=NROW(hist))

	min <- min(start(t))
	offset <- -1*min(start(t))
	pb <- txtProgressBar(min=0, max=NROW(g), style=3)
	for (i in 1:length(g)) {
		setTxtProgressBar(pb, i)
		# access by name
		#hist[as.character(start(t[i,]):end(t[i,]))] <-  hist[as.character(start(t[i,]):end(t[i,]))] + 1
		#histMat[i,(start(t[i,])+offset+1):(end(t[i,])+offset+1)] <- 1
		#histMat[i,(start(t[i,])+offset+1):(end(t[i,])+offset+1)] <- 1
		hist[(start(t[i,])+offset+1):(end(t[i,])+offset+1)] <- hist[(start(t[i,])+offset+1):(end(t[i,])+offset+1)] + 1
	}
	#hist <- colSums(histMat)
	names(hist) <- min(start(t)):max(end(t))

	if (!is.na(filename))
		png(filename)

	if (scale) {
		ylim <- c(0, length(t))
		nT <- length(t)

		plot(min(start(t)):max(end(t)), hist, type="l", xlim=c(-60000,60000), ylim=ylim, xlab="Relative to TSS", ylab="Frequency")
		abline(v=30000, col="blue", lty=2)
		total <- 60000*nT
		block <- 30000*nT
		FP5prime <- trap.rule(1:30000, hist[as.character(-29999:0)])
		TP <- trap.rule(1:30000, hist[as.character(0:29999)])
		PostTTS <- trap.rule(1:30000, hist[as.character(29999:59998)])
		TN <- 30000*nT - FP5prime
		accuracy <- (TP + TN)/ total
		cat("\n")
		print(paste("FP(5'):", round(FP5prime/block,2), "TP:", round(TP/block,2), "PostTTS:",
			round(PostTTS/block,2), "accuracy:", round(accuracy,2)))
		abline(v=-30000, col="blue", lty=2)
	} else {
		plot(min(start(t)):max(end(t)), hist, type="l", ylim=c(0, length(t)))
	}
	abline(v=0, col="blue", lty=2)

	if (!is.na(filename))
	dev.off()

	#title(paste(chrom, strand))
	print(paste("No of transcripts:", length(t)))
	#return(hist)
}

