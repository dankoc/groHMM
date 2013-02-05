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
##	WriteWiggle
##	Date: 2009-07-17
##
##	Writes a wiggle file suitable for uploading to the UCSC genome browser!
##
########################################################################
WriteWiggle <- function(p, file, str="N", size=50, normCounts=1, sep.chrom=FALSE, reverse=FALSE, track.type.line=FALSE, debug=FALSE) { #color="0,0,0", OtherOptions="", 
	F <- WindowAnalysis(p=p, str=str, ssize=size, debug=debug)
	CHR <- as.character(names(F))

	## If we are not separating, prepare the file before the loop
	if(!sep.chrom) {
		filename <- paste(file,".wig", sep="")
		unlink(filename)
		if(track.type.line) {
			write(paste("track type=wiggle_0", sep=""), file=filename, append=TRUE)
		}
	}

	for(i in 1:NROW(CHR)) {
		if(debug) {
			print(paste("Writing:", CHR[i]))
		}

		## If we are separating, prepare the file during the loop!
		if(sep.chrom) {
			filename <- paste(file,CHR[i],".wig", sep="")
			unlink(filename)
			if(track.type.line) {
				write(paste("track type=wiggle_0", sep=""), file=filename, append=TRUE)
			}
		}

		write(paste("fixedStep chrom=", CHR[i], " start=1 step=", size, " span=", size, sep=""), file=filename, append=TRUE)
		if(!reverse) {
			write(as.real(F[[CHR[i]]]*normCounts), file=filename, ncolumns=1, append=TRUE)
		}
		else {
			write(as.real(-1*F[[CHR[i]]]*normCounts), file=filename, ncolumns=1, append=TRUE)
		}
	}

}

