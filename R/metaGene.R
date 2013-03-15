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
##	metaGene
##	Date: 2009-05-27
##
##	Returns a histogram of the number of reads in each section of a
##	moving window centered on a certain feature.
##
##	2009-07-08 -- Name changed to MetaGene.  Formerly, MovingWindow.  
##			Name changed to avoid confusion with new WindowAnalysis function, which is much more general.
##
##	Arguments:
##	features	-> GRanges whose length is 1 
##	(deprecated: f	-> data.frame of: CHR, START, STRAND.)
##	reads	-> GRanges 
##	(deprecated: p	-> data.frame of: CHR, START, END, STRAND.)
##	size	-> The size of the moving window.
##	up	-> Distance upstream of each f to align and histogram.
##	down	-> Distance downstream of each f to align and histogram (NULL).
##
##	Addumptions:
##	(1) 
##
##	TODO: 
##	(1) Add feature accomidating a hemaphrodyte strand (i.e. any ChIP-seq MNase Seq, etc.). 
##	(2) ...
##
########################################################################

#' Returns a histogram of the number of reads in each section of a moving window centered on a certain feature.
#'
#' @param features A GRanges object representing a set of genomic coordinates.  The meta-plot will be centered on the start position.
#' @param reads A GRanges object representing a set of mapped reads.
#' @param size The size of the moving window.
#' @param up Distance upstream of each features to align and histogram Default: 1 kb.
#' @param down Distance downstream of each features to align and histogram Default: same as up.
#' @param debug If set to TRUE, provides additional print options. Default: FALSE
#' @return A vector representing the 'typical' signal centered on a point of interest.
#' @author Charles G. Danko and Minho Chae
metaGene <- function(features, reads, size, up=1000, down=up, debug=FALSE) {
	# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
	features <- features[order(as.character(seqnames(features)), start(features)),]

	if(is.null(down)) {
		down <- up
	}

	C <- sort(unique(as.character(seqnames(features))))
	H <- rep(0,(up + down + 1))
	for(i in 1:NROW(C)) {
		if(debug) {
			print(C[i])
		}

		# Which KG?  prb?
		indxF   <- which(as.character(seqnames(features)) == C[i])
		indxPrb <- which(as.character(seqnames(reads)) == C[i])

		if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
			# Type coersions.
			FeatureStart 	<- start(features[indxF,])
			FeatureStr	<- as.character(strand(features[indxF,]))
			PROBEStart 	<- start(reads[indxPrb,])
			PROBEEnd 	<- end(reads[indxPrb,])
			PROBEStr	<- as.character(strand(reads[indxPrb,]))
			size		<- as.integer(size)
			up			<- as.integer(up)
			down		<- as.integer(down)

			# Set dimensions.
			dim(FeatureStart)	<- c(NROW(FeatureStart), NCOL(FeatureStart))
			dim(FeatureStr)		<- c(NROW(FeatureStr), 	 NCOL(FeatureStr))
			dim(PROBEStart) 	<- c(NROW(PROBEStart), 	 NCOL(PROBEStart))
			dim(PROBEEnd) 		<- c(NROW(PROBEEnd), 	 NCOL(PROBEEnd))
			dim(PROBEStr)		<- c(NROW(PROBEStr), 	 NCOL(PROBEStr))

			if(debug) {
				print(paste(C[i],": Counting reads in specified region.",sep=""))
			}
			Hprime <- .Call("HistogramOfReadsByFeature", FeatureStart, FeatureStr, 
							PROBEStart, PROBEEnd, PROBEStr, 
							size, up, down, PACKAGE = "groHMM")

			H <- H + as.integer(Hprime)
		}
		if(debug) {
			print(paste(C[i],": Done!",sep=""))
		}
	}

	return(H)
}
 
########################################################################
##
##	MetaGeneMatrix
##	Date: 2010-08-27
##
##	Returns a matrix of counts.  Rows represent different streches of DNA.
##	Columns represent positions relative to a certain feature.  Summed together,
##  these should be a meta-gene.
##
##	Arguments:
##	f	-> data.frame of: CHR, START, STRAND.
##	p	-> data.frame of: CHR, START, END, STRAND.
##	size	-> The size of the moving window.
##	up	-> Distance upstream of each f to align and histogram.
##	down	-> Distance downstream of each f to align and histogram (NULL).
##
##	Assumptions: Same as MetaGene
##
########################################################################

#' Returns a matrix, with rows representing read counts across a specified gene, or other features of interest.
#'
#' @param features A GRanges object representing a set of genomic coordinates. 
#' @param reads A GRanges object representing a set of mapped reads. 
#' @param size The size of the moving window.
#' @param up Distance upstream of each f to align and histogram Default: 1 kb.
#' @param down Distance downstream of each f to align and histogram Default: same as up.
#' @param debug If set to TRUE, provides additional print options. Default: FALSE
#' @return A vector representing the 'typical' signal across genes of different length.
#' @author Charles G. Danko and Minho Chae
metaGeneMatrix <- function(features, reads, size= 50, up=1000, down=up, debug=FALSE) {
	# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
	features <- features[order(as.character(seqnames(features)), start(features)),]

	C <- sort(unique(as.character(seqnames(features))))
	H <- NULL
	for(i in 1:NROW(C)) {
		if(debug) {
			print(C[i])
		}

		# Which KG?  prb?
		indxF   <- which(as.character(seqnames(features)) == C[i])
		indxPrb <- which(as.character(seqnames(reads)) == C[i])

		if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
			# Type coersions.
			FeatureStart 	<- start(features[indxF,])
			FeatureStr	<- as.character(strand(features[indxF,]))
			PROBEStart 	<- start(reads[indxPrb,])
			PROBEEnd 	<- end(reads[indxPrb,])
			PROBEStr	<- as.character(strand(reads[indxPrb,]))
			size		<- as.integer(size)
			up			<- as.integer(up)
			down		<- as.integer(down)

			# Set dimensions.
			dim(FeatureStart)	<- c(NROW(FeatureStart), NCOL(FeatureStart))
			dim(FeatureStr)		<- c(NROW(FeatureStr), 	 NCOL(FeatureStr))
			dim(PROBEStart) 	<- c(NROW(PROBEStart), 	 NCOL(PROBEStart))
			dim(PROBEEnd) 		<- c(NROW(PROBEEnd), 	 NCOL(PROBEEnd))
			dim(PROBEStr)		<- c(NROW(PROBEStr), 	 NCOL(PROBEStr))

			if(debug) {
				print(paste(C[i],": Counting reads in specified region.",sep=""))
			}
			Hprime <- .Call("MatrixOfReadsByFeature", FeatureStart, FeatureStr, 
							PROBEStart, PROBEEnd, PROBEStr, 
							size, up, down, PACKAGE = "groHMM")
			H <- rbind(H, Hprime)
		}
		if(debug) {
			print(paste(C[i],": Done!",sep=""))
		}
	}

	return(H)
}



########################################################################
##
##	metaGene_nL
##	Date: 2010-07-05
##
##	Returns a histogram of the number of reads in each section of a
##	moving window of variable size across genes.
##
##	Arguments:
##	f	-> data.frame of: CHR, START, END, STRAND.
##	p	-> data.frame of: CHR, START, END, STRAND.
##	n_windows	-> The resolution of the MetaGene -- i.e. the number of moving windows to break it into..
##
##	Assumptions:
##	(1) Gene list should be ordered!  
##	(2) Gene list should be pretty short, as most of the processing and looping over genes is currently done in R.
##
##	TODO: 
##	(1) Write C function. 
##	(2) ...
##
########################################################################

#' Returns a histogram of the number of reads in each section of a moving window of variable size across genes.
#'
#' @param features A GRanges object representing a set of genomic coordinates. 
#' @param reads A GRanges object representing a set of mapped reads. 
#' @param n_windows The number of windows to break genes into.
#' @param debug If set to TRUE, provides additional print options. Default: FALSE
#' @return A vector representing the 'typical' signal across genes of different length.
#' @author Charles G. Danko and Minho Chae
metaGene_nL <- function(features, reads, n_windows=1000, debug=FALSE) {
	# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
	features <- features[order(as.character(seqnames(features)), start(features)),]

	C <- sort(unique(as.character(seqnames(features))))
	H <- rep(0,n_windows)
	for(i in 1:NROW(C)) {
		if(debug) {
			print(C[i])
		}

		# Which KG?  prb?
		indxF   <- which(as.character(seqnames(features)) == C[i])
		indxPrb <- which(as.character(seqnames(reads)) == C[i])

		if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
			# Type coersions.
			FeatureStart 	<- start(features[indxF,])
			FeatureEnd 	<- end(features[indxF,])
			FeatureStr	<- as.character(strand(features[indxF,]))
			PROBEStart 	<- start(reads[indxPrb,])
			PROBEEnd 	<- end(reads[indxPrb,])
			PROBEStr	<- as.character(strand(reads[indxPrb,]))

			# Set dimensions.
			dim(FeatureStart)	<- c(NROW(FeatureStart), NCOL(FeatureStart))
			dim(FeatureStr)		<- c(NROW(FeatureStr), 	 NCOL(FeatureStr))
			dim(PROBEStart) 	<- c(NROW(PROBEStart), 	 NCOL(PROBEStart))
			dim(PROBEEnd) 		<- c(NROW(PROBEEnd), 	 NCOL(PROBEEnd))
			dim(PROBEStr)		<- c(NROW(PROBEStr), 	 NCOL(PROBEStr))

			for(iFeatures in 1:NROW(FeatureStart)) {
				ws <- (FeatureEnd[iFeatures]-FeatureStart[iFeatures])/n_windows ## This WILL be an interger.
				if(debug) {
					print(paste(C[i],": Counting reads in specified region:",
							FeatureStart[iFeatures],"-",FeatureEnd[iFeatures],sep=""))
					print(paste(C[i],": Window size:",
							FeatureStart[iFeatures],"-",FeatureEnd[iFeatures],sep=""))
					print(paste(C[i],": End-Start:",
							FeatureEnd[iFeatures]-FeatureStart[iFeatures],sep=""))
				}
				DataByOne <- .Call("WindowAnalysis", PROBEStart, PROBEEnd, PROBEStr, FeatureStr[iFeatures],
								as.integer(1), as.integer(1), 
								FeatureStart[iFeatures], FeatureEnd[iFeatures], 
								PACKAGE = "groHMM")

				if(debug) {
					print(paste("DataByOne size:",NROW(DataByOne)))
				}

				## This seems almost immeidate on my pentium M machine.
				Hprime <- unlist(lapply(1:NROW(H), function(i) {
					indx <- ceiling(ws*(i-1)+1):ceiling(ws*i)
					return(sum(DataByOne[indx]))
				}))

				## Reverse it for "-" strands.
				if(FeatureStr[iFeatures] == "-")
					Hprime <- rev(Hprime)

				H <- H + Hprime/ws ## Put in units of: number of reads/base
			}
		}
		if(debug) {
			print(paste(C[i],": Done!",sep=""))
		}
	}

	return(H)
}


########################################################################
##
##	averageGene
##	Date: 2010-12-03
##
##	Returns the average profile of tiling array probe intensity values or wiggle-like count data centered on a set of genomic positions.
##
##	Arguments:
##	Peaks		-> data.frame of: CHR, CENTER, STRAND. (note that STRAND is currenly not supported, and does nothing).
##	ProbeData	-> data.frame of: CHR, CENTER, VALUE
##	bins		-> The bins of the meta gene -- i.e. the number of moving windows to break it into.
##
##	TODO: 
##	(1) Implement support for a Peaks$starnd (
##	(2) ...
##
########################################################################

#' Returns the average profile of tiling array probe intensity values or wiggle-like count data centered on a set of genomic positions (specified by 'Peaks').
#'
#' @param ProbeData Data.frame representing chromosome, window center, and a value.
#' @param Peaks Data.frame representing chromosome, and window center.
#' @param size Numeric.  The size of the moving window. Default: 50 bp.
#' @param bins The bins of the meta gene -- i.e. the number of moving windows to break it into. Default +/- 1kb from center.
#' @return A vector representing the 'typical' signal centered on the peaks of interest.
#' @author Charles G. Danko and Minho Chae
averagePlot <- function(ProbeData, Peaks, size=50, bins= seq(-1000,1000,size)) {

	## For each chromsome.  
	ProbeData$minDist <- rep(999)
	for(chr in unique(Peaks[[1]])) {

	  ##   Get the set of data that fits.  
	  indxPeaks <- which(Peaks[[1]] == chr)
	  
	  ## The '$' ensures that chrom is end of line.  Otherwise, grep for chr1 returns chr10-chr19 as well.
	  ## Should work fine, even when chromsome is simply "chrN".
	  indxAffxProbes <- grep(paste(chr,"$", sep=""), ProbeData[[1]], perl=T)

	  ##   Calculate the minimum distance between the probe and the vector over all features ... 
	  ProbeData$minDist[indxAffxProbes] <- unlist(lapply(indxAffxProbes, function(x) {
			## TODO: For strand to switch it, just multiply by strand here.
			return((ProbeData[x,2] - Peaks[indxPeaks,2])[which.min(abs(ProbeData[x,2] - Peaks[indxPeaks,2]))])}))
	}

	## Make bins.  Take averages and all that...
	means <- unlist(lapply(c(1:NROW(bins)), function(i){mean(ProbeData[(ProbeData$minDist >= (bins[i]-size) & ProbeData$minDist < bins[i]),3])}))
	return(data.frame(windowCenter= bins+(size/2), means))
}


#' runMetaGene Runs Metagene analysis. 
#'
#' Returns a histogram of the number of reads in each section of a moving window centered on a certain feature for both sense and anntisense directions.  It has option for subsampling.
#'
#' @param features GRanges whose length is 1, i.e., Transcription Start Site (TSS) 
#' @param reads GRanges of reads. 
#' @param size Numeric.  The size of the moving window. Default: 100
#' @param up Numeric. Distance upstream of each fgr to align and histogram. Default: 10000
#' @param down Numeric. Distance downstream of each fgr to align and histogram.  If NULL, down is same as up. Default: NULL
#' @param normCounts Numeric.   Reads are multiplied by normCounts.Thhe size of the moving window. Default: 1
#' @param sampling Logical.  If TRUE, subsampling of Metagene is used.  Default: FALSE
#' @param nSampling Numeric. Number of subsampling.  Default: 1000
#' @param debug Logical. If set to TRUE, provides additional print options. Default: FALSE 
#' @return List of vectors representing the 'typical' signal centered on the genomic features of interest.
#' @author Minho Chae
runMetaGene <- function(features, reads, size=100, up=10000, down=NULL, normCounts=1, sampling=FALSE, nSampling=1000, debug=FALSE) {
	if (sampling) {
		plus <- samplingMetaGene(features=features, reads=reads, size=size, up=up, down=down, normCounts=normCounts, 
			nSampling=nSampling, debug=debug)
    	} else {
		plus <- metaGene(features=features, reads=reads, size=size, up=up, down=down, debug=debug)
		plus <- plus / NROW(features)
		plus <- plus*normCounts / NROW(reads)
    	}

	fgrRev <- features
	strand(fgrRev) <- rev(strand(features))
	f <- data.frame(as.character(seqnames(fgrRev)), start(fgrRev), as.character(strand(fgrRev)))
	if (sampling) {
		minus <- samplingMetaGene(features=features, reads=reads, size=size, up=up, down=down, normCounts=normCounts, 
			nSampling=nSampling, debug=debug)
	} else {
		minus <- metaGene(features=features, reads=reads, size=size, up=up, down=down, debug=debug)
		minus <- minus / NROW(features)
		minus <- minus*normCounts / NROW(reads)
	}
	return(list(sense=plus, antisense=minus))
}


samplingMetaGene <- function(features, reads, size=100, up=10000, down=NULL, normCounts=1, nSampling=1000, debug=FALSE) {
	if (is.null(down))
		down <- up

	read_vector <- NULL
	tss_vector <- NULL
	samplingSize <- round(NROW(features)*.1)

	print("Preparing...")
	pb <- txtProgressBar(min=0, max=NROW(features), style=3)
	for(i in 1:NROW(features)) {
		tss_vector <- rbind(tss_vector, metaGene(features=features[i,], reads=reads, size=size, up=up, down=down))
		setTxtProgressBar(pb, i)
	}
	cat("\n")

	print("Sampling...")
	pb <- txtProgressBar(min=0, max=nSampling, style=3)
	M <- matrix(0, nrow=nSampling, ncol=(up+down+1))
	for(i in 1:nSampling) {
		setTxtProgressBar(pb, i)
		sampInx <- sample(1:NROW(features), size=samplingSize, replace=TRUE)
		M[i,] <- apply(tss_vector[sampInx,], 2, sum)
	}
	cat("\n")

	result <- apply(M, 2, median)*normCounts / NROW(reads)
	result <- result/samplingSize
	return(result)
}
