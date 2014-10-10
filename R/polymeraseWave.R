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

##
##  TODO: Re-factor to one function, allowing the model version to be specified 
##  as an argument (Charles).
##
##
##

#'  Given GRO-seq data, identifies the location of the polymerase 'wave' in 
#'  up- or down- regulated genes.  
#'
#'  The model is a three state hidden Markov model (HMM).  States represent: 
#'  (1) the 5' end of genes upstream of the transcription start site, 
#'  (2) upregulated sequence, and (3) the 3' end of the gene through the 
#'  polyadenylation site.  
#'
#'  The model computes differences in read counts between the two conditions.
#'  Differences are assumed fit a functional form which can be specified by 
#'  the user (using the emissionDistAssumption argument).  
#'  Currently supported functional forms include a normal distribution 
#' (good for GRO-seq data prepared using the circular ligation protocol), 
#' a gamma distribution (good for 'spikey' ligation based GRO-seq data), 
#' and a long-tailed normal+exponential distribution was implemented, but 
#' never deployed.
#'
#'  Initial parameter estimates are based on initial assumptions of 
#'  transcription rates taken from the literature.  Subsequently all parameters
#'  are fit using Baum-Welch expetation maximization.
#'
#'  Reference: Danko CG, Hah N, Luo X, Martins AL, Core L, Lis JT, Siepel A, 
#'  Kraus WL. Signaling Pathways Differentially Affect RNA Polymerase II 
#'  Initiation, Pausing, and Elongation Rate in Cells. Mol Cell. 
#'  2013 Mar 19. doi:pii: S1097-2765(13)00171-8. 10.1016/j.molcel.2013.02.015. 
#' 
#'  Arguments:
#'  @param reads1 Mapped reads in time point 1.
#'  @param reads2 Mapped reads in time point 2.
#'  @param genes A set of genes in which to search for the wave.
#'  @param approxDist The approximate position of the wave.  
#'  Suggest using 2000 [bp/ min] * time [min], for mammalian data.
#'  @param size The size of the moving window. Suggest using 50 for direct 
#'  ligation data, and 200 for circular ligation data.  Default: 50.
#'  @param upstreamDist The amount of upstream sequence to include 
#'  Default: 10 kb.
#'  @param TSmooth Optimonally, outlying windows are set a maximum value 
#'  over the inter-quantile interval, specified by TSmooth.  
#'  Reasonable value: 20; Default: NA (for no smoothing).  Users are encouraged
#'  to use this parameter ONLY in combination with the normal distribution 
#'  assumptions.
#'  @param NonMap Optionally, un-mappable positions are trated as missing data.
#'  NonMap passes in the list() structure for un-mappable regions.
#'  @param prefix Optionally, writes out png images of each gene examined for 
#'  a wave.  'Prefix' denotes the file prefix for image names written to disk.
#'  Users are encouraged to create a new directory and write in a full path.
#'  @param emissionDistAssumption Takes values "norm", "normExp", and "gamma".
#'  Specifies the functional form of the 'emission' distribution for states 
#'  I and II (i.e. 5' of the gene, and inside of the wave).  
#'  In our experience, "gamma" works best for highly-variable 'spikey' data, 
#'  and "norm" works for smooth data.  As a general rule of thumb, "gamma" 
#'  is used for libraries made using the direct ligation method, and "norm" 
#'  for circular ligation data.  Default: "gamma".
#'  @param finterWindowSize Method returns 'quality' information for 
#'  each gene to which a wave was fit.  Included in these metrics are several 
#'  that define a moving window.  The moving window size is specified by 
#'  filterWindowSize.  Default: 10 kb.
#'  @param limitPCRDups If true, counts only 1 read at each position 
#'  with >= 1 read.  NOT recommended to set this to TRUE.  Defulat: FALSE.
#'  @param returnVal Takes value "simple" (default) or "alldata". "simple" 
#'  returns a data.frame with Pol II wave end positions.  "alldata" returns all
#'  of the availiable data from each gene, including the full posterior 
#'  distribution of the model after EM.
#'  @param debug If TRUE, prints error messages.
#'  @return Returns either a data.frame with Pol II wave end positions, 
#'  or a List() structure with additional data, as specified by returnVal.
#'  @author Charles G. Danko
#'  @examples
#'  genes <- GRanges("chr7", IRanges(2394474,2420377), strand="+", 
#'  SYMBOL="CYP2W1", ID="54905") 
#'  reads1 <- as(readGAlignments(system.file("extdata", "S0mR1.bam",
#'                              package="groHMM")), "GRanges")
#'  reads2 <- as(readGAlignments(system.file("extdata", "S40mR1.bam",
#'                              package="groHMM")), "GRanges")
#'  approxDist <- 2000*10
#'  # Not run:
#'  # pw <- polymeraseWave(reads1, reads2, genes, approxDist)
##  Given GRO-seq data, identifies the location of the polymerase wave in up- 
##  or down-regulated genes.  This version is based on a full Baum-Welch EM 
##  implementation.
##
##  This is a three state HMM -- initial state representing the intergenic 
##  region 5' of a gene, the second representing the initially upregulated 
##  region, and the third representing the remaining sequence of a gene. 
##
##  We assume that upstream region is intergenic, and thus its emmission 
##  distriubtion is assumed to be a constant, set based on a representative 
##  intergenic region.  This is accomidated in my [1,*) framework by keeping 
##  the vairence constant, and scaling the mean for each gene.
##
## Test with GREB1:     chr2:11,591,693-11,700,363
## GREB1 <- data.frame(chr="chr2", start=11591693, end=11700363, str="+")
polymeraseWave <- function(reads1, reads2, genes, approxDist, size=50, 
    upstreamDist= 10000, TSmooth=NA, NonMap=NULL, prefix=NULL, 
    emissionDistAssumption= "gamma", finterWindowSize=10000, 
    limitPCRDups=FALSE, returnVal="simple", debug=TRUE) {
    if(debug) {
        message("Analyzing windows")
    }   

    genes <- as.data.frame(genes)
    genes <- genes[,c("seqnames", "start", "end", "strand", "SYMBOL", "ID")]

    Fp1 <- windowAnalysis(reads=reads1, strand="+", windowSize=size, 
        limitPCRDups=limitPCRDups)
    Fp2 <- windowAnalysis(reads=reads2, strand="+", windowSize=size, 
        limitPCRDups=limitPCRDups)
    Fm1 <- windowAnalysis(reads=reads1, strand="-", windowSize=size, 
        limitPCRDups=limitPCRDups)
    Fm2 <- windowAnalysis(reads=reads2, strand="-", windowSize=size, 
        limitPCRDups=limitPCRDups)
    sizeP1 <- NROW(reads1)
    sizeP2 <- NROW(reads2)
    expCounts <- mean(NROW(reads1),NROW(reads2))

    ANS <- rep(-1, NROW(genes))
    ENDwave <- rep(-1,NROW(genes))
    STRTwave <- rep(-1,NROW(genes))
    KLdivFinal <- rep(-1,NROW(genes))
    KLdivPar <- rep(-1,NROW(genes))
    minWindLTMed <- rep(FALSE,NROW(genes))
    minMeanWindLTMed <- rep(FALSE,NROW(genes))
    nstates<-as.integer(3)    # number of states in HMM.
    unmap <- NA
    
    ## Possible return value.
    dataList <- list()

## Run the model separately on each gene.
    for(i in seq_along(NROW(genes))) {
        geneData <- list()

        if(debug) {
            message("Starting HMM: ", genes[i,5])
        }

    ############################################################################
    #### Define the gene in terms of the windowed size.
        ## Pull the data for the gene.
        if(genes[i,4] == "+") {
            start <- floor((genes[i,2]-upstreamDist)/size)
            end   <- ceiling(genes[i,3]/size)
            emis1  <- 
                (as.numeric(Fp1[[ as.character(genes[i,1]) ]]))[c(start:end)]
                #/sizeP1*expCounts
            emis2  <- 
                (as.numeric(Fp2[[ as.character(genes[i,1]) ]]))[c(start:end)]
                #/sizeP2*expCounts
        }
        else {
            start <- floor(genes[i,2]/size)
            end   <- ceiling((genes[i,3]+upstreamDist)/size)
            emis1  <- 
                rev((as.integer(Fm1[[ as.character(genes[i,1]) ]]))
                    [c(start:end)])#/sizeP1*expCounts
            emis2  <- 
                rev((as.integer(Fm2[[ as.character(genes[i,1]) ]]))
                    [c(start:end)])#/sizeP2*expCounts
        }
    
        ## Scale to a minimum of 1 read at each position (for fitting Gamma). 
        gene  <- as.numeric(emis1 - emis2)
        if(emissionDistAssumption == "gamma") { 
            ## Leave centered on 0 for the norm_exp/norm emission functions
          gene  <- gene +(-1)*(min(gene))+1 
          ## Must translate points if gamma distributed (gamma undefined <0).
        }
        
        if(is.double(TSmooth)) { ## Interperts it as a fold over the inter 
                                 ## quantile interval to filter.
            message("TSmooth is.integer:", TSmooth)
            medGene <- median(gene)
            iqrGene <- IQR(gene)
            gene[(medGene-gene)>(TSmooth*(iqrGene+1))] <- 
                medGene-(TSmooth*(iqrGene+1))
            gene[(gene-medGene)>(TSmooth*(iqrGene+1))] <- 
                medGene+(TSmooth*(iqrGene+1))
        } else if(!is.na(TSmooth)) {
           gene  <- smooth(gene, kind=TSmooth)
        }

        #write.table(gene, "TMP.gene.Rflat")

#       uTrans<- as.integer(ceiling((upstreamDist)/size))
##      Make the initial guess +5kb --> approxDist.
        uTrans<- as.integer(ceiling((upstreamDist-5000)/size))
        iTrans<- as.integer(ceiling((upstreamDist+approxDist)/size))

        ## Run Baum-Welch
        if(debug) message("initial guess:", uTrans, iTrans, NROW(gene))
        counter <- 0

    ###########################################################################
    ## Calculate a moving average and moving max.

        # For each point.  Left of point is defined as P, right is defined as Q.
        # From min(gene):max(gene).
        # Calculate histogram of points
#      if(!is.null(prefix)) { ## Slow, but required for MinOfMax/Avg filter(s).
        MovMeanSpd <- finterWindowSize#10000

        KLdiv <- rep(0,NROW(gene))
        KS    <- rep(0,NROW(gene))
        Means <- rep(0,NROW(gene))

        MovMean  <- rep(0,NROW(gene))
        MovMax   <- rep(0,NROW(gene))
        dMovMean <- rep(0,NROW(gene))

        for(k in c(2:(NROW(gene)-2))) {
            left  <- gene[c(1:k)]
            right <- gene[c((k+1):NROW(gene))]

            LeftHist  <- hist(left, breaks=c((min(gene, na.rm=TRUE)-1)
                :(max(gene, na.rm=TRUE)+1)), plot=FALSE)
            RightHist <- hist(right, breaks=c((min(gene, na.rm=TRUE)-1)
                :(max(gene, na.rm=TRUE)+1)), plot=FALSE)

            minD <- 0.000001
            KLdiv[k]  <- sum((LeftHist$density*log(
                        (LeftHist$density+minD)/(RightHist$density+minD))))
            KS[k]     <- ks.test(left,right)$statistic[[1]]
            Means[k]  <- mean(left, na.rm=TRUE)-mean(right, na.rm=TRUE)

            MovMean[k] <- mean(gene[max((k-(MovMeanSpd/size)),1)
                :min((k+(MovMeanSpd/size)),NROW(gene))], na.rm=TRUE)
            MovMax[k]  <-  max(gene[max((k-(MovMeanSpd/size)),1)
                :min((k+(MovMeanSpd/size)),NROW(gene))], na.rm=TRUE)
            dMovMean[k]   <- MovMean[k-1] - MovMean[k]
        }
#        }

    ############################################################################
    #### Set up initial paremeter estimates.

        ## Fit transition and initial probabilities.
        tProb  <- as.list(data.frame(
            log(c((1-(1/uTrans)),(1/uTrans),0)),
            log(c(0,(1-(1/(iTrans-uTrans))),(1/(iTrans-uTrans)))), 
            log(c(0, 0, 1))))  # Trans. prob.
        iProb  <- as.double(log(c(1, 0, 0))) # iProb.

        ## Fit initial distribution paremeters for emission probabilities.
        parInt  <- Rnorm(gene[c(1:uTrans)])
        if(is.na(parInt$var) | parInt$var == 0) parInt$var = 0.00001 
        ## Check that the varience of the intergenic state is NOT 0.

        if(emissionDistAssumption == "norm") {
            ePrDist <- c("norm", "norm", "norm") 
            #      ePrDist <- c("norm", "normexp", "normexp")
            parPsi  <- Rnorm(gene[c((uTrans+1):iTrans)])
            #Rnorm.exp(gene[c((uTrans+1):iTrans)], tol=1e-4) #
            parBas  <- Rnorm(gene[c((iTrans+1):NROW(gene))])
            #Rnorm.exp(gene[c((iTrans+1):NROW(gene))], tol=1e-4) #
            ePrVars <- data.frame(c(parInt$mean, 
                sqrt(parInt$var), -1, -1), 
                c(parPsi$mean, sqrt(parPsi$var), -1, -1), 
                c(parBas$mean, sqrt(parBas$var), -1, -1))
        }
        else if(emissionDistAssumption == "normExp") {
            ePrDist <- c("norm", "normexp", "normexp")
            parPsi  <- Rnorm.exp(gene[c((uTrans+1):iTrans)], tol=1e-4) #
            parBas  <- Rnorm.exp(gene[c((iTrans+1):NROW(gene))], tol=1e-4) #
            ePrVars <- data.frame(c(parInt$mean, sqrt(parInt$var), -1, -1),
                    parPsi$parameters, parBas$parameters)
        }
        else if(emissionDistAssumption == "gamma") {
            ePrDist <- c("norm", "gamma", "gamma")
            parPsi  <- RgammaMLE(gene[c((uTrans+1):iTrans)])
            parBas  <- RgammaMLE(gene[c((iTrans+1):NROW(gene))])
            ePrVars <- data.frame(c(parInt$mean, sqrt(parInt$var), -1), 
                c(parPsi$shape, parPsi$scale, -1), 
                c(parBas$shape, parBas$scale, -1))
        }
        else {
          message("emissionDistAssumption should be set to: 'norm', 
            'normExp', or 'gamma'.")
          stopifnot(FALSE) ## Stop here.
        }

#           print(data.frame(c(iMean, parInt$var, -1),c(shape, scale, -1),
#           c(meanBas, stdeBas, -1)))
    ###########################################################################
    ## Finally, choose non-mappable windows and treat them as missing values.
        if(!is.null(NonMap)) {
            windowStarts <- seq(start, end, size) 
            ## Start and end are like chromStart and chromEnd -- 
            ## strand invariant.
            windowEnds <- windowStarts+size
            windowsToSurvey <- data.frame(chrom= rep(as.character(genes[i,1]),
                NROW(windowStarts)), chromStart= as.integer(windowStarts), 
                chromEnd= as.integer(windowEnds), strand= 
                    rep("+",NROW(windowStarts)))
            print(head(windowsToSurvey))
            print(tail(windowsToSurvey))
            unmap <- countMappableReadsInInterval(windowsToSurvey, NonMap)
            if(genes[i,4] == "+") {
                unmap <- rev(unmap)
            }
            gene[unmap/size < 0.10] <- NaN 
            ## Missing data if more than 90% of the windows are unmappable.
    }

    ############################################################################
    ## Now run the HMM.
        if(debug) {
            print(ePrVars)
            print(tProb)
        }

        g <- list()
        g[[1]] <- gene
        ans <- list()
        ans <- tryCatch(.Call("RBaumWelchEM", nstates, g, as.integer(1), 
                ePrDist, ePrVars, tProb, iProb, 0.01, c(TRUE,TRUE,TRUE), 
                c(TRUE, TRUE, TRUE), as.integer(10), TRUE, PACKAGE="groHMM"), 
                error=function(e) e)
                                        ##  Update emis...
        if(NROW(ans) < 3) { 
            ## An error will have a length of 2 (is this guaranteed?!).  
            print("ERROR CAUGHT ON THE C SIDE")
            print(ans) 
            ### Will be a previous iteration of ans...?! Generated a new 'ans'
            ans[[1]] <- NA
            ans[[2]] <- NA
            ans[[3]] <- c(rep(0, NROW(gene)-2), 1, 2)
            ans[[4]] <- NA
            ans[[5]] <- NA
        }
        
        ansVitervi <- ans[[3]][[1]]
        DTs <- max(which(ansVitervi==0))
        DTe <- max(which(ansVitervi==1))

        ### Calculate quality filters... Wrap up.
        KLdivHMM <- 0
        if(debug) print(paste("EM Converged to: ",DTs,DTe,NROW(gene)))
        if((DTs >= 1) & (DTe > 1) & (DTs < DTe) & 
            (DTe < NROW(gene)) & (DTs < NROW(gene))) { 
            ## iff convergest to something useful.
            ## Refit and calculate KL-divergence at that point.
            pINew <- Rnorm(gene[c(1:DTs)])
            
            if(emissionDistAssumption == "norm") {
                pPNew <- Rnorm(gene[c((uTrans+1):iTrans)]) 
                # Rnorm.exp(gene[c((uTrans+1):iTrans)], tol=1e-4) 
                #Rnorm(gene[c((uTrans+1):iTrans)])
                pBNew <- Rnorm(gene[c((DTe+1):NROW(gene))])     
                #Rnorm.exp(gene[c((DTe+1):NROW(gene))], tol=1e-4) #

                ## Estimate KL-divergence.
                PSI <- dnorm(c(min(gene):max(gene)), pBNew$mean, 
                    sqrt(pBNew$var))
                BAS <- dnorm(c(min(gene):max(gene)), pBNew$mean, 
                    sqrt(pBNew$var))
            }
            else if(emissionDistAssumption == "normExp") {
                pPNew <- Rnorm.exp(gene[c((uTrans+1):iTrans)], tol=1e-4) 
                #Rnorm(gene[c((uTrans+1):iTrans)])
                pBNew <- Rnorm.exp(gene[c((DTe+1):NROW(gene))], tol=1e-4) 

                PSI <- (pPNew$parameters[1])*dnorm(c(min(gene):max(gene)), 
                        pPNew$parameters[2], pPNew$parameters[3])+
                        (1-pPNew$parameters[1])*dexp(c(min(gene):max(gene)), 
                        1/pPNew$parameters[4]) ## 1/rate
                BAS <- (pBNew$parameters[1])*dnorm(c(min(gene):max(gene)), 
                        pBNew$parameters[2], pBNew$parameters[3])+
                        (1-pBNew$parameters[1])*dexp(c(min(gene):max(gene)), 
                        1/pBNew$parameters[4]) ## 1/rate
            }
            else if(emissionDistAssumption == "gamma") {
                ## Refit and calculate KL-divergence at that point.
                pPNew <- RgammaMLE(gene[c((DTs+1):DTe)])
                pBNew <- RgammaMLE(gene[c((DTe+1):NROW(gene))])

                ## Estimate KL-divergence.
                PSI <- dgamma(c(min(gene):max(gene)), shape=pPNew$shape, 
                    scale=pPNew$scale)
                BAS <- dgamma(c(min(gene):max(gene)), pBNew$shape, pBNew$scale)
            }
            minD2 <- 1e-300
            KLdivHMM <- sum((PSI*log((PSI+minD2)/(BAS+minD2))))

            ANS[i] <- (DTe-DTs)*size
            STRTwave[i] <- DTs*size
            ENDwave[i] <- DTe*size
            KLdivFinal[i] <- KLdiv[DTe] 
            ##KLdivHMM## Try switching to the empirical KL-div.
            KLdivPar[i] <- KLdivHMM

            ## Calculates min/max and min/avg filters.
            medDns <- median(MovMax[c(max((which(ansVitervi == 1))+
                    round(MovMeanSpd/size)):NROW(MovMax))])
            minMax <- min(MovMax[c(min(which(ansVitervi == 1))
                    :max(which(ansVitervi == 1)))])
            minWindLTMed[i] <- (medDns < minMax) 
            ## True if min(wave) > med(wave.upstream)
            avgDns <- median(MovMean[c(max((which(ansVitervi == 1))+
                    round(MovMeanSpd/size)):NROW(MovMean))])
            minAvg <- min(MovMean[c(min(which(ansVitervi == 1))
                    :max(which(ansVitervi == 1)))])
            minMeanWindLTMed[i] <- (avgDns < minAvg) 
            ## True if min(wave) > med(wave.upstream)
            
        ## Construct the return value...
        geneData[[1]] <- gene[c(1:DTs)]
        geneData[[2]] <- gene[c((DTs+1):DTe)]
        geneData[[3]] <- gene[c((DTe+1):NROW(gene))]
        geneData[[4]] <- emis1 ## Value of the gene... c1.
        geneData[[5]] <- emis2 ## Value of the gene... c2.
        geneData[[6]] <- ans[[4]] ## Matrix of posteriors.
        geneData[[7]] <- ans[[5]] 
        ## Posteriors of a transition from state 2->3.        
        geneData[[8]] <- unmap 
        dataList[[i]] <- geneData
        }
    }
    
    returnDF <- data.frame(StartWave= STRTwave, EndWave= ENDwave, Rate= ANS, 
        KLdiv= KLdivFinal, KLdivParametric= KLdivPar, minOfMax= minWindLTMed, 
        minOfAvg= minMeanWindLTMed, ID= genes[,5], ExternalID= genes[,6])

    if(returnVal == "simple") {
     return(returnDF)
    }
    if(returnVal == "alldata") {
     dataList[[NROW(genes)+1]] <- returnDF
     return(dataList)
    }
} 
