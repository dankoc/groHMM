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
##	MetaGene
##	Date: 2009-05-27
##
##	Returns a histogram of the number of reads in each section of a
##	moving window centered on a certain feature.
##
##	2009-07-08 -- Name changed to MetaGene.  Formerly, MovingWindow.  
##			Name changed to avoid confusion with new WindowAnalysis function, which is much more general.
##
##	Arguments:
##	f	-> data.frame of: CHR, START, STRAND.
##	p	-> data.frame of: CHR, START, END, STRAND.
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
MetaGene <- function(f, p, size, up, down=NULL, debug=FALSE) {
	if(is.null(down)) {
		down <- up
	}

	C <- sort(as.character(unique(f[[1]])))
	H <- rep(0,(up + down + 1))
	for(i in 1:NROW(C)) {
		if(debug) {
			print(C[i])
		}

		# Which KG?  prb?
		indxF   <- which(as.character(f[[1]]) == C[i])
		indxPrb <- which(as.character(p[[1]]) == C[i])

		if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
			# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
			ord <- order(f[indxF,2])

			# Type coersions.
			FeatureStart 	<- as.integer(f[indxF,2][ord])
			FeatureStr	<- as.character(f[indxF,3][ord])
			PROBEStart 	<- as.integer(p[indxPrb,2])
			PROBEEnd 	<- as.integer(p[indxPrb,3])
			PROBEStr	<- as.character(p[indxPrb,4])
			size		<- as.integer(size)
			up		<- as.integer(up)
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
							size, up, down, PACKAGE = "GROseq")

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
MetaGeneMatrix <- function(f, p, size, up, down=NULL, debug=FALSE) {
	if(is.null(down)) {
		down <- up
	}

	C <- sort(as.character(unique(f[[1]])))
	H <- NULL
	for(i in 1:NROW(C)) {
		if(debug) {
			print(C[i])
		}

		# Which KG?  prb?
		indxF   <- which(as.character(f[[1]]) == C[i])
		indxPrb <- which(as.character(p[[1]]) == C[i])

		if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
			# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
			ord <- order(f[indxF,2])

			# Type coersions.
			FeatureStart 	<- as.integer(f[indxF,2][ord])
			FeatureStr	<- as.character(f[indxF,3][ord])
			PROBEStart 	<- as.integer(p[indxPrb,2])
			PROBEEnd 	<- as.integer(p[indxPrb,3])
			PROBEStr	<- as.character(p[indxPrb,4])
			size		<- as.integer(size)
			up		<- as.integer(up)
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
							size, up, down, PACKAGE = "GROseq")

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
##	MetaGene_nL
##	Date: 2010-07-05
##
##	Returns a histogram of the number of reads in each section of a
##	moving window of variable size across genes.
##
##	Arguments:
##	f	-> data.frame of: CHR, START, END, STRAND.
##	p	-> data.frame of: CHR, START, END, STRAND.
##	res	-> The resolution of the MetaGene -- i.e. the number of moving windows to break it into..
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
MetaGene_nL <- function(f, p, res=1000, debug=FALSE) {
	C <- sort(as.character(unique(f[[1]])))
	H <- rep(0,res)
	for(i in 1:NROW(C)) {
		if(debug) {
			print(C[i])
		}

		# Which KG?  prb?
		indxF   <- which(as.character(f[[1]]) == C[i])
		indxPrb <- which(as.character(p[[1]]) == C[i])

		if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
			# Order -- Make sure, b/c this is one of our main assumptions.  Otherwise violated for DBTSS.
			ord <- order(f[indxF,2])

			# Type coersions.
			FeatureStart 	<- as.integer(f[indxF,2][ord])
			FeatureEnd 	<- as.integer(f[indxF,3][ord])
			FeatureStr	<- as.character(f[indxF,4][ord])
			PROBEStart 	<- as.integer(p[indxPrb,2])
			PROBEEnd 	<- as.integer(p[indxPrb,3])
			PROBEStr	<- as.character(p[indxPrb,4])

			# Set dimensions.
			dim(FeatureStart)	<- c(NROW(FeatureStart), NCOL(FeatureStart))
			dim(FeatureStr)		<- c(NROW(FeatureStr), 	 NCOL(FeatureStr))
			dim(PROBEStart) 	<- c(NROW(PROBEStart), 	 NCOL(PROBEStart))
			dim(PROBEEnd) 		<- c(NROW(PROBEEnd), 	 NCOL(PROBEEnd))
			dim(PROBEStr)		<- c(NROW(PROBEStr), 	 NCOL(PROBEStr))

			for(iFeatures in 1:NROW(FeatureStart)) {
				ws <- (FeatureEnd[iFeatures]-FeatureStart[iFeatures])/res ## This WILL be an interger.
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
								PACKAGE = "GROseq")

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
##	AverageGene
##	Date: 2010-12-03
##
##	Returns the average profile of tiling array probe intensity values or wiggle-like count data centered on a set of genomic positions.
##
##	Arguments:
##	Peaks		-> data.frame of: CHR, CENTER, STRAND. (note that STRAND is currenly not supported, and does nothing).
##	ProbeData	-> data.frame of: CHR, CENTER, VALUE
##	bins		-> The bins of the meta gene -- i.e. the number of moving windows to break it into..
##
##	TODO: 
##	(1) Implement support for a Peaks$starnd (
##	(2) ...
##
########################################################################
AveragePlot <- function(ProbeData, Peaks, size=50, bins= seq(-1000,1000,size)) {

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
