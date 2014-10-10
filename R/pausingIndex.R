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

################################################################################
##
##  4/24/2012
##
##  Thanks to Andre Martins for the following two functions, useful for 
##  calculating the confidence interval for a ratio of Poisson random variables.
## 
##  Based on: Ederer F, Mantel N (1974); Confidence Limits on the Ratio of Two 
## Poisson Variables.  AMERICAN JOURNAL OP EPIDEMIOLOGY.
##
############################################################################### 
approx.ratio.CI <- function(x1, x2, alpha=0.05) {
  t = qnorm(1 - alpha/2)
  n = x1 + x2
  zp = (t^2/2 + x1 + 1/2)^2 - ((n + t^2)/n) * (x1 + 1/2)^2
  zn = (t^2/2 + x1 - 1/2)^2 - ((n + t^2)/n) * (x1 - 1/2)^2
  
  a = (t^2/2 + x1 + 1/2 + sqrt(zp)) / (n + t^2/2 - x1 - 1/2 - sqrt(zp))
  b = (t^2/2 + x1 - 1/2 - sqrt(zn)) / (n + t^2/2 - x1 + 1/2 + sqrt(zn))

  return(c(b, a))
}

approx.ratios.CI <- function(num.counts, denom.counts, alpha=0.05) {
  stopifnot(length(num.counts) == length(denom.counts))
  N = length(num.counts)

  result = matrix(data=0, nrow=N, ncol=2)

  for (i in 1:N)
    result[i, ] = approx.ratio.CI(num.counts[i], denom.counts[i], alpha)

  return(result)
}

#' Returns the pausing index for different genes.  TODO: DESCRIBE THE PAUSING 
#' INDEX.
#'
#' Supports parallel processing using mclapply in the 'parallel' package.  
#' To change the number of processors, use the argument 'mc.cores'.
#'
#'  @param features A GRanges object representing a set of genomic coordinates.
#'  @param reads A GRanges object representing a set of mapped reads.
#'  @param size The size of the moving window.
#'  @param up Distance upstream of each f to align and histogram.
#'  @param down Distance downstream of each f to align and histogram (NULL).
#'  @param UnMAQ Data structure representing the coordinates of all un-mappable
#'  regions in the genome.
#'  @param debug If set to TRUE, provides additional print options. 
#'  Default: FALSE
#'  @param ... Extra argument passed to mclapply
#'  @return Returns a data.frame of the pausing indices for the input genes.
#'  @author Charles G. Danko and Minho Chae.
#'  @return Returns the pausing index for different genes.
#'  @examples
#'  features <- GRanges("chr7", IRanges(2394474,2420377), strand="+") 
#'  reads <- as(readGAlignments(system.file("extdata", "S0mR1.bam", 
#'              package="groHMM")), "GRanges")
#'  ## Not run:
#'  # pi <- pausingIndex(features, reads)
##
##  Arguments:
##  f   -> data.frame of: CHR, START, END, STRAND.
##  p   -> data.frame of: CHR, START, END, STRAND.
##  size    -> The size of the moving window.
##  up  -> Distance upstream of each f to align and histogram.
##  down    -> Distance downstream of each f to align and histogram (NULL).
##  UnMAQ   -> Vector of integers representing the coordinates of all 
##  un-MAQable regions in the genome.
##
##  Assumptions:
##  (1) 
##
##  TODO: 
##  (1) Write C function. 
##  (2) Implement: ...
##      (a)  promWind-> As an alternative to up/down, specify a specific window
##      in which to check for the pause-peak.
##      Uses the same format as "f", expect that features are mapped with 
##      respect to the start of f, where 0 indicates the start of 
##      transcription, and negative numbers specify upstream sequence.
##      This is likely useful for rate; perhaps for identifying internal 
##      paused-peaks...
pausingIndex <- function(features, reads, size=50, up=1000, down=1000, 
    UnMAQ=NULL, debug=FALSE, ...) {
    # make sure reads are sorted
    reads <- reads[order(as.character(seqnames(reads)), start(reads)), ]
    f <- data.frame(chrom=as.character(seqnames(features)), 
        start=as.integer(start(features)), end=as.integer(end(features)), 
        strand=as.character(strand(features)))
    if ("symbol" %in% names(mcols(features))){
        f <- cbind(f, symbol=features$symbol) 
    } else {
        f <- cbind(f, symbol=GeneID <- as.character(seq_len(NROW(f))))
    }
    p <- data.frame(chrom=as.character(seqnames(reads)), 
        start=as.integer(start(reads)), end=as.integer(end(reads)), 
        strand=as.character(strand(reads)))

    C <- sort(as.character(unique(f[[1]])))
    Pause <- rep(0,NROW(f))
    Body  <- rep(0,NROW(f))
    Fish  <- rep(0,NROW(f))
    GeneID <- rep("", NROW(f))
    CIl  <- rep(0,NROW(f))
    CIu  <- rep(0,NROW(f))

    ## Pass back information for the fisher test...
    PauseCounts <- rep(0, NROW(f))
    BodyCounts  <- rep(0, NROW(f))
    UpCounts    <- rep(0, NROW(f))
    UgCounts    <- rep(0, NROW(f))

    size        <- as.integer(size)
    up          <- as.integer(up)
    down        <- as.integer(down)

    ###### Calculate PLUS and MINUS index, for DRY compliance.
    PLUS_INDX <- which(f[[4]] == "+")
    MINU_INDX <- which(f[[4]] == "-")

    ###### Identify TSS -- Start for '+' strand, End for '-' strand.
    if(debug) {
        message("Calculating TSS and gene ends for each gene based 
            on strand information.")
    }
    c_tss_indx <- rep(0,NROW(f))
    c_tss_indx[PLUS_INDX] <- 2
    c_tss_indx[MINU_INDX] <- 3
    c_tss <- unlist(lapply(c(1:NROW(f)), function(x) { f[x, c_tss_indx[x]] }))

    ###### Now calculate left and right position for gene body, based 
    ### on '+' or '-'.
    ### Calculate gene end.  Gene start is contiguous with the coordinates 
    ###for the promoter.
    c_gene_end_indx <- rep(0,NROW(f))
    c_gene_end_indx[PLUS_INDX] <- 3
    c_gene_end_indx[MINU_INDX] <- 2
    c_gene_end <- unlist(lapply(c(1:NROW(f)), function(x) { 
                f[x,c_gene_end_indx[x]] }))

    ### Assign left and right.
    gLEFT   <- rep(0,NROW(c_tss))
    gRIGHT  <- rep(0,NROW(c_tss))
    gLEFT[PLUS_INDX]    <- c_tss[PLUS_INDX] + down
    gRIGHT[PLUS_INDX]   <- c_gene_end[PLUS_INDX]
    gLEFT[MINU_INDX]    <- c_gene_end[MINU_INDX]
    gRIGHT[MINU_INDX]   <- c_tss[MINU_INDX] - down

    ## Run parallel version.
    mcp <- mclapply(c(1:NROW(C)), pausingIndex_foreachChrom, C=C, f=f, p=p, 
                    gLEFT=gLEFT, gRIGHT=gRIGHT, c_tss=c_tss, 
                    size=size, up=up, down=down, UnMAQ=UnMAQ, debug=debug, ...)
   
    ## Unlist and re-order values for printing in a nice data.frame.
    for(i in 1:NROW(C)) {
        # Which KG?  prb?
        indxF   <- which(as.character(f[[1]]) == C[i])
                indxPrb <- which(as.character(p[[1]]) == C[i])

        if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
            Pause[indxF][mcp[[i]][["ord"]]]     <- mcp[[i]][["Pause"]]
            Body[indxF][mcp[[i]][["ord"]]]  <- mcp[[i]][["Body"]]
            Fish[indxF][mcp[[i]][["ord"]]]  <- mcp[[i]][["Fish"]]
            GeneID[indxF][mcp[[i]][["ord"]]]    <- mcp[[i]][["GeneID"]]
            
            PauseCounts[indxF][mcp[[i]][["ord"]]] <- mcp[[i]][["PauseCounts"]]
            BodyCounts[indxF][mcp[[i]][["ord"]]]  <- mcp[[i]][["BodyCounts"]]
            UpCounts[indxF][mcp[[i]][["ord"]]]    <- mcp[[i]][["UpCounts"]]
            UgCounts[indxF][mcp[[i]][["ord"]]]    <- mcp[[i]][["UgCounts"]]

            CIl[indxF][mcp[[i]][["ord"]]] <- mcp[[i]][["CIl"]]
            CIu[indxF][mcp[[i]][["ord"]]] <- mcp[[i]][["CIu"]]
        }
    }

    return(data.frame(Pause= Pause, Body= Body, Fisher= Fish, GeneID= GeneID, 
        CIlower=CIl, CIupper=CIu, PauseCounts= PauseCounts, 
        BodyCounts= BodyCounts, uPCounts= UpCounts, uGCounts= UgCounts))
}

pausingIndex_foreachChrom <- function(i, C, f, p, gLEFT, gRIGHT, c_tss, size, 
    up, down, UnMAQ, debug) {
        if(debug) {
            message(C[i])
        }

        # Which KG?  prb?
        indxF   <- which(as.character(f[[1]]) == C[i])
        indxPrb <- which(as.character(p[[1]]) == C[i])

        if((NROW(indxF) >0) & (NROW(indxPrb) >0)) {
            # Order -- Make sure, b/c this is one of our main assumptions.  
            # Otherwise violated for DBTSS.
            Ford <- order(f[indxF,2])
            Pord <- order(p[indxPrb,2])

            # Type coersions.
            FeatureStart    <- as.integer(gLEFT[indxF][Ford])
            FeatureEnd  <- as.integer(gRIGHT[indxF][Ford])
            FeatureStr  <- as.character(f[indxF,4][Ford])
            FeatureTSS  <- as.integer(c_tss[indxF][Ford])

            PROBEStart  <- as.integer(p[indxPrb,2][Pord])
            PROBEEnd    <- as.integer(p[indxPrb,3][Pord])
            PROBEStr    <- as.character(p[indxPrb,4][Pord])

            # Set dimensions.
            dim(FeatureStart)   <- c(NROW(FeatureStart), NCOL(FeatureStart))
            dim(FeatureEnd)     <- c(NROW(FeatureEnd),   NCOL(FeatureEnd))
            dim(FeatureTSS)     <- c(NROW(FeatureTSS),   NCOL(FeatureTSS))
            dim(FeatureStr)     <- c(NROW(FeatureStr),   NCOL(FeatureStr))
            dim(PROBEStart)     <- c(NROW(PROBEStart),   NCOL(PROBEStart))
            dim(PROBEEnd)       <- c(NROW(PROBEEnd),     NCOL(PROBEEnd))
            dim(PROBEStr)       <- c(NROW(PROBEStr),     NCOL(PROBEStr))

        ## Run the calculations on the gene.
            ## Calculate the maximal 50 bp window.
            if(debug) {
                message(C[i],": Counting reads in pause peak.")
            }
            HPause <- .Call("NumberOfReadsInMaximalSlidingWindow", 
                        FeatureTSS, FeatureStr, PROBEStart, PROBEEnd, 
                        PROBEStr, size, up, down, PACKAGE = "groHMM")

            ## Run the calculate on the gene body...
            if(debug) {
                message(C[i],": Counting reads in gene.")
            }
            HGeneBody <- .Call("CountReadsInFeatures", FeatureStart, 
                FeatureEnd, FeatureStr, PROBEStart, PROBEEnd, PROBEStr, 
                PACKAGE = "groHMM")

        ## Get size of gene body
            Difference <- FeatureEnd-FeatureStart
            Difference[Difference < 0] <- 0 
            ## Genes < 1kb, there is no suitable area in the body of the gene.

        ## Calculate UN-MAQable regions...
            if(!is.null(UnMAQ)) {

                ## Count start index.
                chr_indx <- which(UnMAQ[[1]][[1]] == C[i])
                CHRSIZE <- as.integer(UnMAQ[[1]][[2]][chr_indx])
                CHRSTART <- as.integer(0)
                if(chr_indx > 1) {  ## Running on 1:0 gives c(1, 0)
                    CHRSTART <- as.integer( 
                        sum(UnMAQ[[1]][[2]][
                            c(1:(chr_indx-1))
                        ]) +1)
                }

                if(debug) {
                    message(C[i],": Counting unMAQable regions.")
                    message("CHRSIZE:", CHRSIZE, "CHRSTART:", CHRSTART)
                }
                
                ## Count unMAQable regions, and size of everything ... 
                nonmappable <- .Call("CountUnMAQableReads", 
                    FeatureStart, FeatureEnd, as.integer(UnMAQ[[2]]), 
                    CHRSTART, CHRSIZE, PACKAGE = "groHMM")

                ## Adjust size of gene body.
                Difference <- Difference - nonmappable + 1 
                ## Otherwise, get -1 for some.

                if(debug) {
                    print(head(nonmappable))
                    print(as.integer(head(Difference)))
                }
            }

        ## Now use Fisher's Exact.
            if(debug) {
                message(C[i],": Using Fisher's exact.")
            }
            # Make uniform reads.
            Up <- round(HPause + HGeneBody)*(size)/(size+Difference) 
            ## Expted in pause == 
            ##((total reads)/ (total size) [reads/ base]) * 
            ## size [reads/ pause window]
            Ug <- round(HPause + HGeneBody)*(Difference)/(size+Difference) 
            ## Expted reads in body == 
            ## ((total reads)/ (total size) [reads/ base]) * 
            ## (gene size) [reads/ gene body]
            HFish <- unlist(lapply(c(1:NROW(Ford)), function(x) {
               fisher.test(
                  data.frame(
                c(HPause[x], round(Up[x])), 
                c(HGeneBody[x], round(Ug[x]))
                  )
               )$p.value
            } ))

        ## Make return values.
            Pause_c     <- as.double(HPause/size)
            Body_c  <- as.double((HGeneBody+1)/Difference) 
            ## 6-5-2012 -- Add a pseudocount here, 
            ## forcing at least 1 read in the gene body.   
            Fish_c  <- as.double(HFish)
            GeneID_c    <- as.character(f[indxF,5][Ford])
            
            PauseCounts_c <- HPause
            BodyCounts_c  <- HGeneBody
            UpCounts_c    <- round(Up)
            UgCounts_c    <- round(Ug)

            aCI <- approx.ratios.CI(HPause, HGeneBody)
            scaleFactor <- Difference/size ## Body/ pause, must be 1/ PI units.
            CIl_c <- as.double(aCI[,1]*scaleFactor)
            CIu_c <- as.double(aCI[,2]*scaleFactor)
            
            return(list(Pause= Pause_c, Body= Body_c, Fish= Fish_c, 
                GeneID= GeneID_c, PauseCounts= PauseCounts_c, 
                BodyCounts= BodyCounts_c, UpCounts= UpCounts_c, 
                UgCounts= UgCounts_c, CIl= CIl_c, CIu= CIu_c, ord= Ford))
        }
    return(integer(0))
}
