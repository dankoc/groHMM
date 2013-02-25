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

#` writeWiggle writes a wiggle track suitable for uploading to the UCSC genome browser.
#` Currently only supports writing a fixed-step wiggle.
#`
#` @param reads GenomicRanges object representing the position of reads mapping in the genome.
#` @param file Specifies the file prefix for the output wiggle.
#` @param strand Takes values of "+", "-", or "N".  Computes Writes a wiggle on the speicified strand.  "N" denotes collapsing reads on both strands.  Default: "N".
#` @param size Size of non-overlapping windows to write. Default: 50 bp.
#` @param normCounts A normalization factor correcting for library size or other effects.  For example, total mappible read counts might be a reasonable value.  Default: 1 (i.e. no normalization).
#` @param sep.chrom If set to TRUE, will write a separate wiggle file for each chromosome.  Default: FALSE
#` @param reverse If set to TRUE, multiplies values by -1.  Used for reversing GRO-seq data on the negative (-) strand. Default: FALSE
#` @param track.type.line If set to TRUE, prints a header identifying the file as a wiggle.  Necessary to upload a custom track to the UCSC genome browser.  Default: TRUE
#` @param debug If set to TRUE, provides additional print options. Default: FALSE
#` @return Writes a wiggle file to the specified file.
writeWiggle <- function(reads, file, strand="N", size=50, normCounts=1, sep.chrom=FALSE, reverse=FALSE, track.type.line=TRUE, debug=FALSE) { #color="0,0,0", OtherOptions="", 

	## Error checking. ... 
	if(!(strand=="N"|strand=="+"|strand=="-")) {
	  stop("Strand should be specified as '+', '-', or 'N'.")
	}
	
	F <- windowAnalysis(reads=reads, strand=strand, ssize=size, debug=debug)
	CHR <- as.character(names(F))

	## If we are not separating, prepare the file before the loop
	if(!sep.chrom) {
		filename <- paste(file,".wig", sep="")
		unlink(filename)
		if(track.type.line) {
			write(paste("track type=wiggle_0", sep=""), file=filename, append=TRUE)
		}
	}

    cat("Writing...\n")
    pb <- txtProgressBar(min=0, max=NROW(CHR), style=3)
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
        setTxtProgressBar(pb, i)
	}
    cat("\n")
}

