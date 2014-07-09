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
 *  Source code written for GROseq package by Charles Danko.
 *
 *  2009-05-07 Started this file, writing 
 *
 ******************************************************************************/

/**************************************************************
 *
 *  Associates a vector of genomic featuers (e.g. genes, CpG islands, etc.) 
 *  with a table of sequence reads.
 *
 **************************************************************/
#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// HMM algorithm entry points.
SEXP Rviterbi(SEXP emi, /*SEXP seqLength,*/ SEXP nEmis, SEXP nstates, 
    SEXP emiprobDist, SEXP emiprobVars, SEXP tprob, SEXP iprob);

SEXP RBaumWelchEM(/*SEXP nSeq,*/ SEXP nstates, SEXP emi, 
    /*SEXP seqLength*/ SEXP nEmis, SEXP emiprobDist, SEXP emiprobVars, 
    SEXP tprob, SEXP iprob, SEXP threshold, SEXP updatetrans, SEXP updateemis, 
    SEXP output, SEXP verbose);

SEXP HistogramOfReadsByFeature(SEXP FeatureStart, SEXP FeatureStrand, 
                SEXP ReadStart, SEXP ReadEnd, SEXP ReadStrand, 
                SEXP size, SEXP up, SEXP down);

SEXP CountUnMAQableReads(SEXP FeatureStart, SEXP FeatureEnd, SEXP UnMAQ, 
    SEXP offset, SEXP sizeofchr);

SEXP WindowAnalysis(SEXP ProbeStart, SEXP ProbeEnd, SEXP ProbeStrand, 
    SEXP CheckStrand, SEXP windowsize, SEXP stepsize, SEXP startposition, 
    SEXP endposition);

SEXP DecayAlgorithm(SEXP COUNTS, SEXP DECAY);
SEXP getTranscriptPositions(SEXP Transform, SEXP Threshold, SEXP WindowSize);
SEXP vect2bed(SEXP Transform, SEXP WindowSize);

SEXP CountReadsInFeatures(SEXP Feature_Start, SEXP Feature_End, 
    SEXP Feature_Chr, SEXP ProbeStart, SEXP ProbeEnd, SEXP ProbeChr);

SEXP AssociateRegionWithFeatures(SEXP Feature_Start, SEXP Feature_End, 
    SEXP ProbeStart, SEXP ProbeLength);

// In MLEfit
SEXP RgammaMLE(SEXP n, SEXP sumxi, SEXP sumlogxi);

/**************************************************************
 *
 *  Register entry points...
 *
 **************************************************************/
void R_init_groHMM(DllInfo *info) {
     R_CallMethodDef callMethods[]  = {
       {"AssociateRegionWithFeatures", 
        (DL_FUNC)&AssociateRegionWithFeatures, 4},
       {"CountReadsInFeatures", (DL_FUNC)&CountReadsInFeatures, 6},
       {"Rviterbi", (DL_FUNC)&Rviterbi, 7},
       {"RBaumWelchEM", (DL_FUNC)&RBaumWelchEM, 12},
       {"HistogramOfReadsByFeature", (DL_FUNC)&HistogramOfReadsByFeature, 8},
       {"CountUnMAQableReads", (DL_FUNC)&CountUnMAQableReads, 5},
       {"WindowAnalysis", (DL_FUNC)&WindowAnalysis, 8},
       {"DecayAlgorithm", (DL_FUNC)&DecayAlgorithm, 2},
       {"getTranscriptPositions", (DL_FUNC)&getTranscriptPositions, 3},
       {"vect2bed", (DL_FUNC)&vect2bed, 2},
       {"RgammaMLE", (DL_FUNC)&RgammaMLE, 3},
       {NULL, NULL, 0}
     };

    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

