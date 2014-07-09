/***************************************************************************
**
**   Copyright 2009, 2010, 2011 Charles Danko.
**
**   This program is part of the groHMM R package
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
 *  Source code written for GRO-seq package by Charles Danko.
 *
 *  2009-07-07 More funcitons added for all useful windowing analysis.
 *  2009-05-27 Started this file, specifically used to identify pausing indices 
 *
 ******************************************************************************/

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
 * DecayAlgorithm -- Returns the number of reads in each bin ...
 *
 * 2009-07-08: Wrote this wrapper ...
 *  Assume that DECAY is expressed in units of decay per window.
 *
 * 
 ******************************************************************************/
SEXP DecayAlgorithm(SEXP COUNTS, SEXP DECAY) {

    int *counts = INTEGER(COUNTS);
    double decay = REAL(DECAY)[0];

    // Get the dimensions.
    SEXP DIM1;
    DIM1 = getAttrib(COUNTS,R_DimSymbol);
    int size = INTEGER(DIM1)[0];

    // Construct return values.
    SEXP ModCounts;
    PROTECT(ModCounts = allocVector(REALSXP,size));
    double *modcounts = REAL(ModCounts);

    // Construct return values. Just reutrn the origional vector of counts?!?!
    modcounts[0] = counts[0];
    for(int i=0;i<size-1;i++) {
//      modcounts[i+1] = 
//      ((double)modcounts[i]*decay<counts[i+1])?counts[i+1]:
//          (double)modcounts[i]*decay;
        modcounts[i+1] = (double)modcounts[i]*decay+counts[i+1];
        /*template <class T> const T& max ( const T& a, const T& b ) {
        return (b<a)?a:b;     // or: return comp(b,a)?a:b; for the comp version
        }*/
    }

    UNPROTECT(1);
    return(ModCounts);
}

/******************************************************************************
 *
 * DEPRICATED!!!!!!!!!!!!!!!!!!!!!!!
 * 
 * getTranscriptPositions -- Converts vector of transcript positions to a 
 * bed-like file.
 *
 * DEPRICATED!!!!!!!!!!!!!!!!!!!!!!!
 *
 ******************************************************************************/
SEXP getTranscriptPositions(SEXP Transform, SEXP Threshold, SEXP WindowSize) {

    const int false= 0;
    // const int true= 1; unused

    double *transform = REAL(Transform);
    double threshold = REAL(Threshold)[0];
    int windowsize = INTEGER(WindowSize)[0];

    // Get the dimensions.
    int size = Rf_nrows(Transform);

    // Run through once counting the size ...
    int currentlyInRegion=false;
    int NumberOfRegions=0;
    for(int i=0;i<size;i++) {
        if( (transform[i] >= threshold) && !currentlyInRegion ) {
            currentlyInRegion = !currentlyInRegion;
            NumberOfRegions++;
        }
        else if( (transform[i] < threshold) && currentlyInRegion ) {
            currentlyInRegion = !currentlyInRegion;
        }
    }

    // Construct return values.
    SEXP Regions, RegionStarts, RegionEnds, COL_Names;
    PROTECT(Regions = allocVector(VECSXP, 2));
        PROTECT(COL_Names = NEW_CHARACTER(2));

    SET_VECTOR_ELT(Regions, 0, RegionStarts=allocVector(INTSXP, 
            NumberOfRegions));
    SET_STRING_ELT(COL_Names, 0, mkChar("Start"));

    SET_VECTOR_ELT(Regions, 1, RegionEnds=allocVector(INTSXP, 
            NumberOfRegions));
    SET_STRING_ELT(COL_Names, 1, mkChar("End"));

        setAttrib(Regions, R_NamesSymbol, COL_Names);

    int *starts = INTEGER(RegionStarts);
    int *ends = INTEGER(RegionEnds);

    // Find the starts/ends.
    currentlyInRegion=false;
    int RegionNumber=0;
    // Init to -1.

    // use: (size-1)? -- don't include that last individual window at the end 
    // of the chrom?!?!
    for(int i=0;i<(size);i++) {
        if( (transform[i] >= threshold) && !currentlyInRegion ) {
            currentlyInRegion = !currentlyInRegion;
            if(RegionNumber < NumberOfRegions) {
                starts[RegionNumber] = i*windowsize;
                ends[RegionNumber] = i*windowsize;
            }
            else {
                break;
            }
        }
        else if( (transform[i] < threshold) && currentlyInRegion ) {
            currentlyInRegion = !currentlyInRegion;
            if(RegionNumber < NumberOfRegions) {
                ends[RegionNumber++] = i*windowsize+windowsize;
            }
            else  {
                break;
            }
        }
    }

    UNPROTECT(2);
    return(Regions);
}


/*******************************************************************************
 *
 * vect2bed -- Converts vector of transcript positions to a bed-like file.
 *
 * Will replace getTranscriptPositions
 * 
 ******************************************************************************/
SEXP vect2bed(SEXP Transform, SEXP WindowSize) {

    double *transform = REAL(Transform);
    double currValue=transform[0];
    int windowsize = INTEGER(WindowSize)[0];

    // Get the dimensions.
    int size = Rf_nrows(Transform);

    // Run through once counting the size ...
    int NumberOfRegions=1;  // Include the first region ... write any time it 
                            // switches...
    for(int i=0;i<size;i++) {
        if( transform[i] != currValue ) {
            currValue=transform[i];
            NumberOfRegions++;
        }
    }

    // Construct return values.
    SEXP Regions, RegionStarts, RegionEnds, StateID, COL_Names;
    PROTECT(Regions = allocVector(VECSXP, 3));
    PROTECT(COL_Names = NEW_CHARACTER(3));

    SET_VECTOR_ELT(Regions, 0, RegionStarts=allocVector(INTSXP, 
            NumberOfRegions));
    SET_STRING_ELT(COL_Names, 0, mkChar("Start"));

    SET_VECTOR_ELT(Regions, 1, RegionEnds=allocVector(INTSXP, 
            NumberOfRegions));
    SET_STRING_ELT(COL_Names, 1, mkChar("End"));
    
    SET_VECTOR_ELT(Regions, 2, StateID=allocVector(INTSXP, NumberOfRegions));
    SET_STRING_ELT(COL_Names, 2, mkChar("State"));

    setAttrib(Regions, R_NamesSymbol, COL_Names);

    int *starts = INTEGER(RegionStarts);
    int *ends = INTEGER(RegionEnds);
    int *state = INTEGER(StateID);

    // Find the starts/ends.
    int RegionNumber=0;
    currValue=transform[0];
    starts[RegionNumber]=0*windowsize;
    state[RegionNumber]=transform[0];
    // Init to -1.

    // use: (size-1)? -- 
    // don't include that last individual window at the end of the chrom?!?!
    for(int i=0;i<size;i++) {
        if( transform[i] != currValue ) {
            if(RegionNumber < NumberOfRegions) {
             ends[RegionNumber++] = i*windowsize+windowsize;
             starts[RegionNumber] = i*windowsize;
             state[RegionNumber] = transform[i];
             currValue = transform[i];
            }
            else {
                Rprintf("WARNING! Size of variable EXCEEDED! \
                    It's really a MAJOR PROBLEM!");
                break;
            }
        }
    }

    UNPROTECT(2);
    return(Regions);
}

