########################################################################
##
##      makeConsensusAnnotations
##      Date: 2012-11-27 
## 
##	Makes a consensus annotation which is a non overlapping annotation where multiply covered isoforms are used
##	to represent a genomic ranges for the gene.   This also reduces annotations in a way to save starting site of
##	a gene for overlapping gene annotations.
##
##      Arguments:
##      anno     -> GRanges of annotation
## 	    keytype  -> Keytype to collapse gene annotation, i.e., gene_id 
##
##      Assumptions:
##      (1) 
##
##      TODO: 
##      (1) 
##
########################################################################
#' makeConsensusAnnotations Makes a consensus annotation 
#'
#' Makes a consensus annotation which is a non overlapping annotation where multiply covered isoforms are used
#' to represent a genomic ranges for the gene.   This also reduces annotations in a way to save starting site of
#' a gene for overlapping gene annotations.
#'
#' @param annotations GRanges of annotations to be collapsed. 
#' @param keytype Character.  Keytype to collapse isoforms.  Default: "gene_id"
#' @return GenomicRanges object of reads. 
#' @author Charles G. Danko and Minho Chae
makeConsensusAnnotations <- function(annotations, keytype="gene_id") {
	# First reduce by keytype, collapsing isoforms
	annoByKeytype <- reduceAnnotationByKeytype(annotations, keytype=keytype)
	return(reduceOverlappingAnnotation(annoByKeytype))
}


reduceAnnotationByKeytype <- function(anno, keytype="gene_id") {
	if (!(keytype %in% colnames(mcols(anno))))
		stop("Please supply a keytype argument.")

	# No isoforms 
	annoList <- GenomicRanges::split(anno, GenomicRanges::unlist(mcols(anno)[,keytype]))
	uInx <- elementLengths(annoList) == 1  

	uniGenes <- GenomicRanges::unlist(annoList[names(which(uInx))])
	print(paste("Genes with no isoforms:", NROW(uniGenes)))

	# process genes with isoforms
	isoList <- annoList[!uInx]
	isoGenes <- GenomicRanges::unlist(endoapply(isoList, function(x) x[1,]))
	print(paste("Genes with isoforms:", NROW(isoGenes)))

	print(paste("Reducing by ", keytype, "...", sep=""))
	pb <- txtProgressBar(min=0, max=NROW(isoList), style=3)

	for (i in 1:NROW(isoList)) {
	x <- isoList[[i]] # Make sure they are all same strand
	geneID <- GenomicRanges::unlist(elementMetadata(x)$gene_id)[1]
	nPlus <- sum(strand(x)=="+")
	nMinus <- sum(strand(x)=="-")
	if (nPlus*nMinus != 0) {  # mixed strand
		if ((nPlus > nMinus) || (nPlus == nMinus)) {
			x <- x[strand(x)=="+",]
		} else {
			x <- x[strand(x)=="-",]
		}
	}

	d <- disjoin(x)
	if (NROW(d) == 1) {   # there is only one disjoint interval
		result <- d
	} else {
		o <- findOverlaps(d, x)
		multi <- queryHits(o)[duplicated(queryHits(o))]  # disjoint interval with multiple overlapping (coverage)
		if (NROW(multi) == 0) {
			result <- d[1,]                              # take just any disjoint interval
		} else {
			single <- setdiff(1:NROW(d), multi)              # disjoint interval with single coverage
			if (NROW(single) == 0 ) {
				reduced <- reduce(d)
			} else {
				reduced <- reduce(d[-single,])
			}

			if (NROW(reduced) == 1) {
				result <- reduced
			} else {    # merge reduced by filling gaps
					result <- GRanges(seqnames=seqnames(reduced[1,]),
						ranges=IRanges(start=min(start(reduced)), end=max(end(reduced))),
						strand=strand(reduced[1,]))
			}
		}
	}

	seqnames(isoGenes[geneID]) <- seqnames(result)
	ranges(isoGenes[geneID]) <- ranges(result)
	strand(isoGenes[geneID]) <- strand(result)
	setTxtProgressBar(pb, i)
	}
	cat("\n")

	final <- c(uniGenes, isoGenes)
	return(final[order(as.character(seqnames(final)), start(final))])
}

 # There are still overlapping annotations
reduceOverlappingAnnotation <- function(anno) {
	ol <- findOverlaps(anno, ignoreSelf=TRUE, ignoreRedundant=TRUE)

	intx <- pintersect(anno[queryHits(ol),], anno[subjectHits(ol),])
	oQual <- data.frame( query=queryHits(ol), subject=subjectHits(ol),
			   qOverlap=width(intx)/width(anno[queryHits(ol),]),
			    sOverlap=width(intx)/width(anno[subjectHits(ol),]),
			    overBases=width(intx),
			    similarity=width(intx)/pmax(width(anno[queryHits(ol),]), width(anno[subjectHits(ol),])))

	print(paste("No of total overlappings:", NROW(oQual)))
	# Remove cases where one annotation is totally contained by another, this will be taken care by "reduce" later
	oQual <- oQual[!(oQual$qOverlap == 1 | oQual$sOverlap == 1),]
	if (NROW(oQual) == 0 ) {
		consensus <- anno
	} else {
		print(paste("No of partial containment::", NROW(oQual)))
		cat("Reducing overlapping annotations...\n")
		pb <- txtProgressBar(min=0, max=NROW(oQual), style=3)
		brAnno <- GRanges()

		for (i in 1:NROW(oQual)) {
			overlapped <- anno[as.integer(oQual[i,c("query","subject")]),]
			overlapped <- overlapped[order(start(overlapped))]
			dis <- disjoin(overlapped)  # three fragments
			if (as.character(strand(dis[1,])) == "+") {  # last two combined
				one <- overlapped[1,]
				ranges(one) <- ranges(dis[1,])
				end(one) <- end(one) - 1
				combined <- overlapped[2,]
				ranges(combined) <- ranges(dis[2,])
				end(combined) <- end(dis[3,])
			} else {        # first two combined
				one <- overlapped[2,]  # use it as a template
				ranges(one) <- ranges(dis[3,])
				start(one) <- start(one) + 1 # Need more gap, otherwised it's going to be reduced
				combined <- overlapped[1,]
				ranges(combined) <- ranges(dis[1,])
				end(combined) <- end(dis[2,])
			}
			brAnno <- c(brAnno, one)
			brAnno <- c(brAnno, combined)
			setTxtProgressBar(pb, i)
		}

        	# Now remove overlapped one from anno and add disjoint ones
        	consensus <- anno[-unique(c(oQual$query, oQual$subject)),]
        	consensus <- c(consensus, brAnno)
        	cat("\n")
	}

	# reduce by saving metadata
	return(reduceBySavingMetadata(consensus))
}

reduceBySavingMetadata <- function(anno) {
	ol <- findOverlaps(anno, ignoreSelf=TRUE, ignoreRedundant=TRUE)
	if (length(ol) == 0)
		return(anno[order(as.character(seqnames(anno)), start(anno))])

	olInx <- sort(unique(c(queryHits(ol), subjectHits(ol))))
	nonOverlapped <- anno[setdiff(1:NROW(anno), olInx),]   # non overlapped annotations
	overlapped <- anno[olInx,]                             # overlapped annotations

	reduced <- reduce(overlapped)

	rol <- findOverlaps(reduced, overlapped)
	rol.df <- data.frame(reduced=queryHits(rol), overlapped=subjectHits(rol))

	curRed <- rol.df[1,"reduced"]
	curMax <- rol.df[1, "overlapped"]
	maxList <- numeric()
	for (i in 2:NROW(rol.df)) {
		newRed <- rol.df[i, "reduced"]
		newOver <- rol.df[i, "overlapped"]
		if (curRed == newRed) {
			if (which.max(c(width(overlapped[curMax,]), width(overlapped[newOver,]))) == 1) {
			# Same 
			} else {
				curMax <- newOver
			}
		} else {  # current reduce index changed
			maxList <- c(maxList, curMax)
			curRed <- newRed
			curMax <- newOver
		}
	}
	# last one
	maxList <- c(maxList, curMax)

	# add metadata
	elementMetadata(reduced) <- elementMetadata(overlapped[maxList,])

	final <- c(nonOverlapped, reduced)
	return(final[order(as.character(seqnames(final)), start(final))])
}


