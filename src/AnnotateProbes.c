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


/*******************************************************************************
 *
 *  Source code written for GROseq package by Charles Danko.
 *
 *  2009-05-07 Started this file, writing 
 *
 ******************************************************************************/

/**************************************************************
 *
 *  Associates a vector of genomic featuers (e.g. genes, CpG islands, etc.) with
 *  a table of sequence reads.
 *
 **************************************************************/
#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/******************************************************************************
 *
 * CountReadsInFeatures -- Count the number of probes in each genomic feature 
 * (i.e. regions).
 *
 *  Arguments:
 *  Feature_Start   -> integer vector representing the start of features along 
 *      the chromosome.
 *  Feature_End -> integer vector representing the end of features along the 
 *      chromosome.
 *  Feature_Chr -> chromosome of each features, either '+' or '-'.
 *  ProbeStart  -> vector of integers that represents the start of each probe 
 *      along the chromosome.
 *  ProbeEnd    -> vector of integers that represents the end of each probe.
 *  ProbeChr    -> chromosome of each read, either '+' or '-'.
 *
 *  Returns:
 *  integer array -- for each Feature, the number of Probes/Reads inside the 
 *      feature, on the same strand.
 *
 *  Assumes:
 *  (1) All are on the same chromosome.
 *  (2) All Features are sorted in ascending order.
 *  (3) All Probes are sorted in ascending order.
 *
 ******************************************************************************/
SEXP CountReadsInFeatures(SEXP Feature_Start, SEXP Feature_End, 
    SEXP Feature_Chr, SEXP ProbeStart, SEXP ProbeEnd, SEXP ProbeChr) {
    int *fSTART = INTEGER(Feature_Start);
    int *fEND = INTEGER(Feature_End);
    int *PS = INTEGER(ProbeStart);
    int *PE = INTEGER(ProbeEnd);

    // Get the dimensions.
    SEXP DIM1, DIM2;
    DIM1 = getAttrib(Feature_Start,R_DimSymbol);
    int NFEATURES = INTEGER(DIM1)[0];
    DIM2 = getAttrib(ProbeStart, R_DimSymbol);
    int NPROBES = INTEGER(DIM2)[0];

    // Construct return values.
    SEXP fID;
    PROTECT(fID = allocVector(INTSXP,NFEATURES));
    int *fcID = INTEGER(fID);

    // Assign probes to a feature
/*  for(int features=0;features<NFEATURES;features++) {
        fcID[features] = 0;
        const char *feature_str = CHAR(STRING_ELT(Feature_Chr, features)); 
        // Should NOT allocate new memory?!
        for(int probes=0;probes<NPROBES;probes++) {
            if((fSTART[features] <= (PE[probes])) && 
                (fEND[features] >= PS[probes]) &&
              (strcmp(feature_str, CHAR(STRING_ELT(ProbeChr, probes))) == 0)) {
                    fcID[features]++;
            }
        }
    }*/

    // Assign probes to a feature.
    // Tested just fine using: RefSeqNH00 <- 
    // CountReadsInInterval(f= Gbed, p= NH00[,c(1:3,6)])
    int counter;
    int prev_counter_start=0; // Start from offset 'o'.
    for(int features=0;features<NFEATURES;features++) {
        fcID[features] = 0;
        const char *feature_str = CHAR(STRING_ELT(Feature_Chr, features)); 
        // Should NOT allocate new memory?!

        // Figure out where to start, w/ some error checking
        if(fSTART[features] > PE[prev_counter_start -1]) 
            counter = prev_counter_start;
        else counter = 0;

        while((fSTART[features] > PE[counter]) && (counter < NPROBES)) 
            counter++; // Find MAQ @ feature start.

        while((  fEND[features] >= PS[counter]) && (counter < NPROBES)) {
            if(strcmp(feature_str, CHAR(STRING_ELT(ProbeChr, counter))) == 0) {
                fcID[features]++; // Detect a decrease.
            }
            prev_counter_start = counter; 
            // Features are in order, so start from here next time.
            counter++;
        }

    }


    UNPROTECT(1);
    return(fID);
}


/*******************************************************************************
*
* Associates probes with genomic features -- or regions.
*
*   Arguments:
*   Feature_Start   -> integer vector representing the start of features along 
*       the chromosome.
*   Feature_End -> integer vector representing the end of features along the 
*       chromosome.
*   ProbeStart  -> vector of integers that represents the start of each probe 
*       along the chromosome.
*   ProbeLength -> vector of probe lengths
*
*   Returns:
*   integer array -- for each ProbeStart, the index of the feature that it is 
*       inside, or NA for none.
*
*   Assumes:
*   The feature table passed as Feature_Start and Feature_End are on the 
*       same chromosome.
*   Feature_Start and Feature_End are of the same length.
*
*   09-05-23 -- Copied this function to GROseq from AffyTiling.
*
*******************************************************************************/
SEXP AssociateRegionWithFeatures(SEXP Feature_Start, SEXP Feature_End, 
    SEXP ProbeStart, SEXP ProbeLength) {
    int *fSTART = INTEGER(Feature_Start);
    int *fEND = INTEGER(Feature_End);
    int *PS = INTEGER(ProbeStart);
    int *PL = INTEGER(ProbeLength);

    // Get the dimensions.
    SEXP DIM1, DIM2;
    DIM1 = getAttrib(Feature_Start,R_DimSymbol);
    int NCPG = INTEGER(DIM1)[0];
    DIM2 = getAttrib(ProbeStart, R_DimSymbol);
    int NPROBES = INTEGER(DIM2)[0];

    // Construct return values.
    SEXP fID;
    PROTECT(fID = allocVector(INTSXP,NPROBES));
    int *fcID = INTEGER(fID);

    // Assign probes to a feature
    for(int prb=0;prb<NPROBES;prb++) {
        fcID[prb] = NCPG+1; // This will add NA when names are added in R.
        for(int i=0;i<NCPG;i++) {
            if((fSTART[i] < (PS[prb] + PL[prb])) && (fEND[i] > PS[prb])) {
                fcID[prb] = (i+1);
                break; // Only 1 feature for each probe!
            }
        }
    }

    UNPROTECT(1);
    return(fID);
}
