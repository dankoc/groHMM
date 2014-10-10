/***************************************************************************
**
**   Copyright 2009, 2010, 2011 Charles Danko.
**
**   This program is part of the GRO-seq R package
**
**   groHMM is free software: you can redistribute it and/or modify it 
**   under the terms of the GNU General Public License as published by 
**   the Free Software Foundation, either version 3 of the License, or  
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful, but 
**   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
**   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
**   for more details.
**
**   You should have received a copy of the GNU General Public License along 
**   with this program.  If not, see <http://www.gnu.org/licenses/>.
**
***************************************************************************/


/******************************************************************************
 *
 *	Source code written for GRO-seq package by Charles Danko.
 *
 *	2009-07-07 More funcitons added for all useful windowing analysis.
 *	2009-05-27 Started this file, specifically used to identify pausing indices 
 *
 ******************************************************************************/


/* using namespace std; */

/* extern "C" { */

#include <R.h> 
#include <S.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*******************************************************************************
 *
 * SlidingWindow -- Count number of reads in a sliding window, along a chomosome
 *  w/ variable step size.
 *
 * 2009-07-07 Written to allow comparison between GRO-seq and POLLII-ChIP-seq --
 *  designed to cover entire chromosomes.
 *
 * Arguments:
 *	(1)  Probe_Start	-> *long int, of size NProbes, indicating the start point 
 *	   of each probe/read.
 *	(2)  Probe_End		-> *long int, of size NProbes, indicating the end point of 
 *	    each probe/read.
 *	(3)  Probe_Strand	-> *char, '+' or '-', of size NProbes, indicating the 
 *	    strand of each probe/read.
 *	(4)  NProbes		-> The number of probe/reads.
 *	(5)  Straind		-> The strand on which to count, pass "N" to ignore 
 *	    strand.
 *	(6)  WindowSize		-> Thesize of the sliding window (bp).
 *	(7)  StepSize		-> The size of the step between each window position (bp).
 *	(8)  StartPosition	-> Chromosome start.
 *	(9)  EndPosition	-> Chromosome end.
 *	(10) InitialIndex	-> Index of the initial point, from which to scan.  
 *	    Use 0 if unsure.
 *		NOTE:  InitialIndex is not yet used, and ignored.
 *
 * Returns:
 *	(1) Vector of counts of reads for each window in the interval.
 *
 * Assumes:
 *	(1) DOES NOT CURRENTLY ASSUME THAT Probes are in order!!
 *	(2) DOES NOT CURRENTLY use probe strand information! 
 *	(3) StepSize > 0.
 *
 ******************************************************************************/

//int *SlidingWindow(	int *Probe_Start, int *Probe_End, SEXP Probe_Strand, 
//    int NProbes, const char *Strand, int WindowSize, int StepSize, 
//    int StartPosition, int EndPosition, int &InitialIndex, int *counts) {

int *SlidingWindow(	int *Probe_Start, int *Probe_End, SEXP Probe_Strand, 
    int NProbes, const char *Strand, int WindowSize, int StepSize, 
    int StartPosition, int EndPosition, int *counts) {

// Make a new *int counts, the proper size.
//	Equal to the number of window START positions in the region.  
//	Depends ONLY on step size -- not on window size.
	int size  = ceil((double)(EndPosition-StartPosition)/(double)StepSize);
	for(int i=0;i<size;i++)
		counts[i] = 0;

// Foreach Probe, (should scale linearly with the number of probes).
	for(int p=0;p<NProbes;p++) {

//	Continue if we are outside the window.
		if(Probe_End[p] < StartPosition) 
			continue;
		if(EndPosition < Probe_Start[p])
			break;


// 	Continue if we are on a different strand AND we are not ignoring strand.
		if( !(strcmp(Strand, "N") == 0) 
			&& !(strcmp(Strand, CHAR(STRING_ELT(Probe_Strand, p))) == 0) )
			continue;

//	Increment each window near each probe.
		// Find w such that probestart is in it.
		// Find probestart ...
		// Proof is in 7/8/09 notebook entry (January 2009 notebook).
		int indxFirst = floor((double)
            (Probe_Start[p]-StartPosition-WindowSize)/(double)StepSize);
		int indxLast  = 
            ceil((double)(Probe_End[p]-StartPosition)/(double)StepSize);

		for(int w=indxFirst;w<=indxLast;w++) {
			if(  	(w >= 0) && (w < size) && 
                // Some error checking!  Next, check both conditions:
				(Probe_End[p] >= (StepSize*w+StartPosition)) && 
				((StepSize*w+WindowSize+StartPosition) >= Probe_Start[p]) ) 
					counts[w]++;
		}
	}

	return counts;
}

/******************************************************************************
 *
 * WindowAnalysis -- Returns the number of reads in each bin ...
 *
 * 2009-07-08: Wrote this wrapper ...
 *
 * 
 ******************************************************************************/
SEXP WindowAnalysis(SEXP ProbeStart, SEXP ProbeEnd, SEXP ProbeStrand, 
    SEXP CheckStrand, SEXP windowsize, SEXP stepsize, SEXP startposition, 
    SEXP endposition) {

	int II = 0;

	int *Probe_Start = INTEGER(ProbeStart);
	int *Probe_End = INTEGER(ProbeEnd);
	int *WindowSize = INTEGER(windowsize);
	int *StepSize = INTEGER(stepsize);
	int *StartPosition = INTEGER(startposition);
	int *EndPosition = INTEGER(endposition);

	// Get the dimensions.
	SEXP DIM1;
	DIM1 = getAttrib(ProbeStart,R_DimSymbol);
	int NProbes = INTEGER(DIM1)[0];

	// Construct return values.
	int size  = 
        ceil((double)(EndPosition[0]-StartPosition[0])/(double)StepSize[0]);
	SEXP COUNTS;
	PROTECT(COUNTS = allocVector(INTSXP,size));
	int *counts = INTEGER(COUNTS);

	// SlidingWindow(	Probe_Start, Probe_End, ProbeStrand, NProbes, 
    //     CHAR(STRING_ELT(CheckStrand, 0)), WindowSize[0], StepSize[0], 
    //    StartPosition[0], EndPosition[0], II, counts);

	SlidingWindow(	Probe_Start, Probe_End, ProbeStrand, NProbes, 
        CHAR(STRING_ELT(CheckStrand, 0)), WindowSize[0], StepSize[0], 
        StartPosition[0], EndPosition[0], counts);

	UNPROTECT(1);
	return(COUNTS);
}


/******************************************************************************
 *
 * MetaSlidingWindow -- Calculate number of reads in a sliding window, given 
 *  genomic coordinates.  Anchor points make the metagene-type analysis 
 *  very convenient!
 *
 * Arguments:
 *	(1) Anchor_Start	-> long int, indicating the start point of the anchor.
 *	(2) Anchor_Strand	-> character, '+' or '-', indicating the strand of the 
 *	    anchor point.
 *	(3) Probe_Start		-> *long int, of size NProbes, indicating the start point 
 *	    of each probe/read.
 *	(4) Probe_End		-> *long int, of size NProbes, indicating the end point of 
 *	    each probe/read.
 *	(5) Probe_Strand	-> *char, '+' or '-', of size NProbes, indicating the 
 *	    strand of each probe/read.
 *	(6) NProbes		-> The number of probe/reads.
 *	(7) WindowSize		-> The size of the sliding window (bp).
 *	(8) mdUpstream*		-> The distance from the sliding point to measure, 
 *	    upstream of the anchor point.
 *	(9) mdDownstream*	-> The distance from the sliding point to measure, 
 *	    downstream of the anchor point.
 *	(10) InitialIndex	-> Index of the initial point, from which to scan.  
 *	    Use 0 if unsure.
 *	(11) ans		-> Pointer to pre-allocated space.  Allocating in this function 
 *	    causes memory problems.
 *
 * Returns:
 *	(1) A vector of counts, each separated by 1bp.  
 *	    Size= mdUpstream+mdDownstream+1;
 *	(2) Return vector [size] is the index, found to be the first index in the 
 *	    round that just finished.
 *		Since probes are assumed to be both orderd and in the same chromosome, 
 *		the index of the next gene is guaranteed to be >= Return[size].  
 *		Meant to be a time-saver.
 *
 * Assumes: 
 *	(1) All positions, etc. are on the same chromosome.  
 *	(2) Positions are all in order!
 *	(3) NO option for stepsize -- stepsize = 1.
 *	(4) This general function takes and returns C-types.
 *
 * Move this function into its own cpp file --- tools, or something like that.
 *
 ******************************************************************************/

int *MetaSlidingWindow(int Anchor_Start, const char *Anchor_Strand, 
			int *Probe_Start, int *Probe_End, SEXP Probe_Strand, int NProbes,
			int WindowSize, int mdUpstream, int mdDownstream, int InitialIndex, 
            int *ans) {

	int First, Last;
	if(Anchor_Strand[0] == '+') {
		First = Anchor_Start - mdUpstream - WindowSize; // Depending on strand!
		Last  = Anchor_Start + mdDownstream + WindowSize;
	}
	else if (Anchor_Strand[0] == '-') {
		First = Anchor_Start - mdDownstream - WindowSize; // Depending on strand!
		Last  = Anchor_Start + mdUpstream + WindowSize;
	}
	else {
		error("Incorrect strand: %s",Anchor_Strand);
	}

	// Find index from InitialIndex, from which to start the search...  
    // Run some error checking.
	int INDX = InitialIndex;
	if((INDX < 0) || (INDX >= NProbes)) INDX = 0; // Error checking ...
	if((INDX > 0) && (Probe_Start[INDX-1] > First)) INDX = 0; 
    // More error checking...

	int size = mdUpstream + mdDownstream + 1; 
    // Length up & down & the zero position (=1).
	for(int i=0;i<(size+1);i++)		// Init all to 0.
		ans[i]=0;

// One loop starting at INDX.  Does not record unless its in the window.
	int indx=0;     // Initialize
	int InWindow=0; // FALSE
	for(int i=INDX;i<NProbes;i++) {
		// If ANYWHERE in the region of interest AND on the same strand as 
        // feature (alternatively, don't count the strand if 'N'); record.
		if((First <= (Probe_End[i])) && (Last >= Probe_Start[i]) &&
		  ((strcmp(Anchor_Strand, CHAR(STRING_ELT(Probe_Strand, i))) == 0) || 
          (strcmp("N", CHAR(STRING_ELT(Probe_Strand, i))) == 0) )) {
			// Put the start of the first INDX.
			if(!InWindow) {
				ans[size] = i;
				InWindow=1; // TRUE
			}

			// Foreach bin that contains a read, Increment!
			for(int j=(Probe_Start[i] - First - 2*WindowSize);
						j<(Probe_End[i] - First -1);j++) {
				if( (j >= 0) && (j < size) ) { // Error check ... 
					// Turn around, if '-' strand, turn around to align 
                    // '+' and '-' together.
					if (Anchor_Strand[0] == '+') indx=j;
					else if (Anchor_Strand[0] == '-') indx=size-j-1;

					ans[indx]++; // RECORD :)
				}
			}
		}
		// Assuming that reads are in order -- 
        // cuts out after it passes the region with reads.
		else if(Probe_Start[i] > Last) {
			return ans;
		}
	}

	return ans; // Only gets here if it hits the end...  
                // In this case, it should be all 0's.
}

/*******************************************************************************
 *
 * CountUnMAQableReads -- Returns a vector of counts over the features...
 *
 * Feature*		-> Vector of information on features start and strand.
 * Read*		-> Vector of information on read start, end, and strand.
 * UnMAQ		-> Int vector representing the genomic coordinates (hg18) of 
 *  non-unique 44mers.
 * offset		-> Sum of the size of all proceeding chromosomes ...
 * sizeofchr		-> Number of un-MAQable regions in the chromosome -- prevents 
 *      reading past the chromosome.
 * 
 *****************************************************************************/
SEXP CountUnMAQableReads(SEXP FeatureStart, SEXP FeatureEnd, SEXP UnMAQ, 
    SEXP offset, SEXP sizeofchr) {

	int *fSTART = INTEGER(FeatureStart);
	int *fEND = INTEGER(FeatureEnd);
	int *MAQ = INTEGER(UnMAQ);
	int *o = INTEGER(offset);
	int *s = INTEGER(sizeofchr);
	int max = o[0]+s[0];

	// Get the dimensions.
	SEXP DIM1;
	DIM1 = getAttrib(FeatureStart,R_DimSymbol);
	int NFEATURES = INTEGER(DIM1)[0];

	// Construct return values.
	SEXP fID;
	PROTECT(fID = allocVector(INTSXP,NFEATURES));
	int *fcID = INTEGER(fID);

	// Assign probes to a feature.  
	int un_counter;
	int prev_un_counter_start=o[0]; // Start from offset 'o'.
	for(int features=0;features<NFEATURES;features++) {
		fcID[features] = 0;

		// Figure out where to start, w/ some error checking
		if(fSTART[features] > MAQ[prev_un_counter_start -1]) un_counter = 
            prev_un_counter_start;
		else un_counter = o[0];

		while((fSTART[features] > MAQ[un_counter]) && (un_counter <= max)) 
			un_counter++; // Find MAQ @ feature start.

		while((  fEND[features] >= MAQ[un_counter]) && (un_counter <= max)) {
			fcID[features]++; // Detect a decrease.
			prev_un_counter_start = un_counter; 
            // Features are in order, so start from here next time.
			un_counter++;
		}

	}

	UNPROTECT(1);
	return(fID);
}


/******************************************************************************
 *
 * HistogramOfReadsByFeature -- Returns a histogram over the feature...
 *
 * Feature.*		-> Vector of information on features start and strand.
 * Read.*		-> Vector of information on read start, end, and strand.
 * size			-> The size of the moving window (bp).
 * up			-> The distance upstream of each feature to generate the 
 *          histogram for.
 * down			-> The distance downstream of each feature to make a histogram 
 *          for.
 *
 ******************************************************************************/

SEXP HistogramOfReadsByFeature(SEXP FeatureStart, SEXP FeatureStrand, 
				SEXP ReadStart, SEXP ReadEnd, SEXP ReadStrand, 
				SEXP size, SEXP up, SEXP down) {

	// Assign variables
	int *gSTART 	= INTEGER(FeatureStart);
	int *PS 	= INTEGER(ReadStart);
	int *PE 	= INTEGER(ReadEnd);

	int SIZE 	= INTEGER(size)[0];
	int UP 		= INTEGER(up)[0];
	int DOWN	= INTEGER(down)[0];

	// Get the dimensions.
	SEXP DIM1, DIM2;
	DIM1 = getAttrib(FeatureStart,R_DimSymbol);
	int NGENES = INTEGER(DIM1)[0];
	DIM2 = getAttrib(ReadStart, R_DimSymbol);
	int NREADS = INTEGER(DIM2)[0];

	// Construct return values.
	int sz = UP + DOWN + 1;
	SEXP counts;
	PROTECT(counts = allocVector(INTSXP,sz));
	int *COUNTS = INTEGER(counts);
	for(int i=0;i<sz;i++) // INIT to 0.
		COUNTS[i] = 0;

	// Foreach gene, run the comparison.
	int *ADD = (int*)R_alloc(sz, sizeof(int));
	int InitIndex = 0;
	for(int i=0;i<NGENES;i++) {
		MetaSlidingWindow(gSTART[i], CHAR(STRING_ELT(FeatureStrand, i)), 
            PS, PE, ReadStrand, NREADS, SIZE, UP, DOWN, InitIndex, ADD);

		// Foreach ADD, sum with 
		for(int j=0;j<sz;j++)
			COUNTS[j] += ADD[j];

		InitIndex = ADD[sz];
	}

	UNPROTECT(1);
	return(counts);
}

/*******************************************************************************
 *
 * MatrixOfReadsByFeature -- Returns a matrix representing counts of reads over 
 *  the feature...
 *
 * Feature.*		-> Vector of information on features start and strand.
 * Read.*		-> Vector of information on read start, end, and strand.
 * size			-> The size of the moving window (bp).
 * up			-> The distance upstream of each feature to generate the 
 *  histogram for.
 * down			-> The distance downstream of each feature to make a histogram 
 *  for.
 *
 ******************************************************************************/

SEXP MatrixOfReadsByFeature(SEXP FeatureStart, SEXP FeatureStrand, 
				SEXP ReadStart, SEXP ReadEnd, SEXP ReadStrand, 
				SEXP size, SEXP up, SEXP down) {

	// Assign variables
	int *gSTART 	= INTEGER(FeatureStart);
	int *PS 	= INTEGER(ReadStart);
	int *PE 	= INTEGER(ReadEnd);

	int SIZE 	= INTEGER(size)[0];
	int UP 		= INTEGER(up)[0];
	int DOWN	= INTEGER(down)[0];

	// Get the dimensions.
	SEXP DIM1, DIM2;
	DIM1 = getAttrib(FeatureStart,R_DimSymbol);
	int NGENES = INTEGER(DIM1)[0];
	DIM2 = getAttrib(ReadStart, R_DimSymbol);
	int NREADS = INTEGER(DIM2)[0];

	// Construct return values.
	int sz = UP + DOWN + 1;
	SEXP counts;
	PROTECT(counts = allocMatrix(INTSXP,NGENES,sz));
	int *COUNTS = INTEGER(counts);
	for(int i=0;i<sz;i++) // INIT to 0.
	  for(int j=0;j<NGENES;j++)
		COUNTS[NGENES*i+j] = 0;

	// Foreach gene, run the comparison.
	int *ADD = (int*)R_alloc(sz, sizeof(int));
	int InitIndex = 0;
	for(int i=0;i<NGENES;i++) {
		MetaSlidingWindow(gSTART[i], CHAR(STRING_ELT(FeatureStrand, i)), 
            PS, PE, ReadStrand, NREADS, SIZE, UP, DOWN, InitIndex, ADD);

		// Foreach ADD, sum with 
		for(int j=0;j<sz;j++)
			COUNTS[NGENES*j+i] += ADD[j];
			
		InitIndex = ADD[sz];
	}

	UNPROTECT(1);
	return(counts);
}

/******************************************************************************
 *
 * NumberOfReadsInMaximalSlidingWindow -- Foreach feature, returns the number 
 *  of reads in the maximal sliding window in range ...
 *
 * Feature.*		-> Vector of information on features start and strand.
 * Read.*		-> Vector of information on read start, end, and strand.
 * size			-> The size of the moving window (bp).
 * up			-> The distance upstream of each feature to generate the 
 *  histogram for.
 * down			-> The distance downstream of each feature to make a histogram 
 *  for.
 *
 ******************************************************************************/

SEXP NumberOfReadsInMaximalSlidingWindow(SEXP FeatureStart, SEXP FeatureStrand, 
				SEXP ReadStart, SEXP ReadEnd, SEXP ReadStrand, 
				SEXP size, SEXP up, SEXP down) {

	// Assign variables
	int *gSTART 	= INTEGER(FeatureStart);
	int *PS 	= INTEGER(ReadStart);
	int *PE 	= INTEGER(ReadEnd);
	int SIZE 	= INTEGER(size)[0];
	int UP 		= INTEGER(up)[0];
	int DOWN	= INTEGER(down)[0];
	SEXP DIM1, DIM2;
	DIM1 = getAttrib(FeatureStart,R_DimSymbol);
	int NGENES = INTEGER(DIM1)[0];
	DIM2 = getAttrib(ReadStart, R_DimSymbol);
	int NREADS = INTEGER(DIM2)[0];

	// Return values.
	SEXP counts;
	PROTECT(counts = allocVector(INTSXP,NGENES));
	int *COUNTS = INTEGER(counts);
	for(int i=0;i<NGENES;i++)
		COUNTS[i] = 0;

	// Foreach gene, run the comparison.
	int sz = UP + DOWN + 1;
	int *ADD = (int*)R_alloc(sz, sizeof(int));

	int InitIndex = 0;
	for(int i=0;i<NGENES;i++) {
		// Run the sliding window.
		MetaSlidingWindow(gSTART[i], CHAR(STRING_ELT(FeatureStrand, i)), 
            PS, PE, ReadStrand, NREADS, SIZE, UP, DOWN, InitIndex, ADD);

		// Foreach ADD, Look for the maximum 
		for(int j=0;j<sz;j++)
			if(COUNTS[i] < ADD[j])
				COUNTS[i] = ADD[j];

		// Apply InitIndex, for speed!
		InitIndex = ADD[sz];
	}

	UNPROTECT(1);
	return(counts);
}
/* } */
