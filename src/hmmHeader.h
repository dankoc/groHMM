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


/*************************************************************
 *
 *  This file is comprised of prototypes and structure 
 *  definitions for the components of the HMM class.
 *
 *  Started: 11-25-2009.
 *
 *************************************************************/

#ifndef HMMHEADER_H
#define HMMHEADER_H

// Assume that: exp(a) + exp(b) == exp(a); if(a-b > APPROX_EXP_VALUE_THRESHOLD)
#include "UsefulValues.h"
#include <R_ext/Applic.h>
#include <R.h>

//  Generalized emission function.
typedef double (**emiss_func)(double value, double* args, int nArgs);
//(double data, double arg1, double arg2, double arg3);

/****************************************
 *
 * Data structures for the various HMM
 * components.
 *
 ****************************************/
 
typedef struct {
  double *log_iProb;    /* (1 x n_states) Log of initial probabilities. */
  double **log_tProb;   
  /* (n_states x n_states) table of transition probabilities. */
  emiss_func log_eProb; 
  /* log emission probability; CHARLES: Changed to an array of unique 
   * functions for each state. */
  double **em_args; 
  /* CHARLES: Add an appropriate container (i.e. **double) to 
   * store vars. for emission prob. */
  int n_states;     /* number of states (excluding initial state). */
  int n_emis; /* number of columns in emissions matrix **data */
} hmm_t;

typedef struct {
  double **forward;
  double **backward;
  double **data;
  hmm_t *hmm;
  int N; /* sequence size; length(data) */
  double log_px;
  double bk_log_px;
} fwbk_t;

// Can we make this more general?  A separate ss_t "member" struct 
// for each distribution?
// Not quite so object-friendly as java :)
// How about separating and using a void type ?!?
typedef struct {
  double sumPiXi;   
  // \sum_{i}(Pi, Xi) --> 
  // (posterior probability of being in state)*(emission value).
  double sumLogPiXi;    
  // \sum_{i}(log(Pi, Xi)) --> 
  // (posterior probability of being in state)*(emission value)..
  double sumPiXiSq; 
  // This is useful in fitting first version of parameters in 
  // update by method of moments.
  double N;     // \sum_{i}(Pi).
} ssGamma;

typedef struct {
  double sumPiXi;   
  // \sum_{i}(Pi, Xi) --> 
  // (posterior probability of being in state)*(emission value).
  double N;     // \sum_{i}(Pi).
} ssPoisson;

typedef struct {
  double sumPiXi;   
  // \sum_{i}(Pi, Xi) --> 
  // (posterior probability of being in state)*(emission value).
  double N;     // \sum_{i}(Pi).
} ssExp;

typedef struct {
  double sumPiXi;
  double sumPiXiSq;
  double N;
} ssNormal;

typedef struct {
  double *ex; // Can just pass into R's cgmin function.
  int containsData;
} ssNormExp;

typedef struct {
  double **totalTransK; 
  // Matrix, indexed: var[states in HMM][training seqeunces]
} ssTransition;

/*********************************************
 * Generalized sufficient statis prototypes.
 *********************************************/

// To be as memory-friendly as possible, 
// updating sufficient stats is done in four pieces.  
// First, an appropriate sufficient_stats struct is allocated.  
// Ideally, a separate struct
//  will be used for each prob. distribution.
typedef void* (**alloc_emis_sstats)(int num); 
// num --> the number of data points.

// Second, update_suffstats_func keeps a running tally of the 
// sufficient stats for a distribution.
typedef void (**update_sstat_func)(int state, int emis_indx, 
    void* ss, fwbk_t fwbk);

// After all sufficient stats are added from all chromosomes, 
// new variables are calculated and applied to the hmm.
typedef void (**update_emiss_func)(int state, void* ss, hmm_t *hmm);

// Free ss_t variables.
typedef void  (**free_emis_sstats)(void* ss);

// Similar paradigm for transition prob.
typedef void* (**alloc_trans_sstats)(int num, int sequences); 
// num --> the number of data points.
typedef void  (**update_trans_SS)(int state, int sequence, void* ss, 
    emiss_func EMI, fwbk_t fwbk);
typedef void  (**update_trans_Prob)(int state, int nSequences, void* ss, 
    hmm_t *hmm);
typedef void  (**free_trans_sstats)(void* ss);

/****************************************
 *
 * Struct containers for EM vars.
 *
 ****************************************/

 typedef struct {
    // Store function pointers for Transition sufficient stats.
    alloc_trans_sstats  AllocTssFunc;
    update_trans_SS    UpdateTssFunc;
    update_trans_Prob   UpdateTPFunc;
    free_trans_sstats    FreeTssFunc;
    void** TransSS;
    
    // Store function pointers for Emission sufficient stats.
    alloc_emis_sstats sstats_alloc;
    update_sstat_func sstats_emis;
    update_emiss_func update_emis;
    free_emis_sstats  free_emis_s;
    void** ss;
    
    // To update, or not to update? That is the question ...
    int *updateTrans;
    int *updateEmis;
 } em_t;

/***************************************
 * 
 * Prototype allocation/free functions.
 *
 ***************************************/
extern double **matrix_alloc(int rows, int cols, int reverse);
extern void matrix_free(double ** matrix, int reverse, int rows);

extern int **imatrix_alloc(int rows, int cols, int reverse);
extern void imatrix_free(int ** matrix, int reverse, int rows);

extern fwbk_t * fwbk_alloc(double **data, int N, hmm_t *hmm); 
// This may be converted into the R entry function?!
extern void fwbk_free(fwbk_t * data);
extern hmm_t *setupHMM(SEXP nstates, SEXP emiprobDist, 
    SEXP emiprobVars, SEXP nEmis, SEXP tprob, SEXP iprob);

/***************************************
 * 
 * Prototype HMM algorithms.
 *
 ***************************************/
extern void viterbi_path(hmm_t hmm, double **data, 
    int seq_len, double **matrix, int **backptr, int *path);
extern void backward(fwbk_t *data);
extern void forward(fwbk_t *data);

/**************
 * 
 *  Emission probability update functions.  
 *
 *  Two functions are used for each supported distribution:
 *   (1) SStats* -- updates the sufficient statistics with the results 
 *      of a new fwbk run.
 *   (2) Update* -- using the sufficient statistics, updates the paremeters 
 *      in hmm.
 *
 **************/
extern void *SSallocNormal(int num);
extern void SStatsNormal(int state, int emis_indx, void* ss, fwbk_t fwbk);
extern void UpdateNormal(int state, void* ss, hmm_t *hmm);
extern void SSfreeNormal(void* ss);

extern void *SSallocNormExp(int num);
extern void SStatsNormExp(int state, int emis_indx, void* ss, fwbk_t fwbk);
extern void UpdateNormExp(int state, void* ss, hmm_t *hmm);
extern void SSfreeNormExp(void* ss);

extern void *SSallocGamma(int num);
extern void SStatsGamma(int state, int emis_indx, void* ss, fwbk_t fwbk);
extern void SStatsGamma_p1(int state, int emis_indx, void* ss, fwbk_t fwbk);
// Used to add 1 to each iteration.
extern void UpdateGamma(int state, void* ss, hmm_t *hmm);
extern void UpdateGamma_SCALE1(int state, void* ss, hmm_t *hmm); 
// Used to fit a constrained gamma where E[mean] = E[var].
// Used to fit a constrained gamma, where E[x] = 1, and shape=1/scale.
extern void UpdateGamma_SHAPEeq1overSCALE(int state, void* ss, hmm_t *hmm);
extern void SSfreeGamma(void* ss);

extern void *SSallocPoisson(int num);
extern void SStatsPoisson(int state, int emis_indx, void* ss, fwbk_t fwbk);
extern void UpdatePoisson(int state, void* ss, hmm_t *hmm);
extern void SSfreePoisson(void* ss);

extern void *SSallocExp(int num);
extern void SStatsExp(int state, int emis_indx, void* ss, fwbk_t fwbk);
extern void UpdateExp(int state, void* ss, hmm_t *hmm);
extern void SSfreeExp(void* ss);

extern void *TransAlloc(int num, int sequences); 
// num --> the number of data points.
extern void  TransUpdate(int state, int sequence, void* ss, 
    emiss_func EMI, fwbk_t fwbk);
extern void  TransUpdateP(int state, int nSequences, void* ss, hmm_t *hmm);
extern void  TransFree(void* ss);

/*************************************
 *
 * Prototype of MLE for Gamma distribution must also appear here... 
 *
 *************************************/
extern int MLEGamma(double N, double SumXi, double SumLogXi, double *shape, 
    double *scale);
extern int MLEGamma_SCALE1(double N, double SumXi, double SumLogXi, 
    double *shape, double *scale);
extern int MLEGamma_SHAPEeq1overSCALE(double N, double SumXi, 
    double SumLogXi, double SumXiSq, double *shape, double *scale);
extern void normal_exp_optimgr(int n, double *par, double *gr, void *ex);
extern double normal_exp_optimfn(int n, double *par, void *ex);

/*************************************************
 *
 * DRY complience functions.
 *
 * Inline functions supposed to be implemented in the header file ?!
 *
 *************************************************/

// Returns the log(prob.) of  being in state=state, 
// at sequence position=position, given the fwbk output= fwbk.
static  double MargainalizeSumLogProbOver(int state, int position, 
    fwbk_t fwbk) {  
    double logPP = 
        fwbk.forward[position][state]+fwbk.backward[position][state]-
            fwbk.log_px;
    return(logPP);
}
  
 /**************
  * expDif -- Robustly returns log(abs(exp(pGr)-exp(pLs))).
  **************/ 
static  double expDif(double pLs, double pGr){
  if(pGr == pLs) 
    return(log(0));
  else if(pGr > pLs)
    return(log(1-exp(pLs-pGr))+pGr);
  else
    return(log(1-exp(pGr-pLs))+pLs);
 }

static  double expSum(double *logValues, int length) {
  double scalefactor, CurrentSum;
  double TotalSum=0;
  scalefactor = logValues[0];
  for(int k=1;k<length;k++)
    scalefactor= max(scalefactor, logValues[k]);

  for(int k=0;k<length;k++) {
    CurrentSum=logValues[k]-scalefactor;
    if(-1*CurrentSum < APPROX_EXP_VALUE_THRESHOLD)
    TotalSum += exp(CurrentSum);
  }

 return(log(TotalSum) + scalefactor);
}

// expSum w/ just two values.  Inelegant, yet saves having to copy 
// two doubles to a new double*.
static  double expSum2(double v1, double v2) {
  double scalefactor = max(v1, v2);
  double retVal= log(exp(v1-scalefactor)+exp(v2-scalefactor))+scalefactor;
  return(retVal);
}

/**************
 * 
 *  Emission probability functions -- using density functions from R 
 *  where possible.  The functions below are set up such that they will 
 *  pass to Rlibrary functions for different distributions.  
 *  Mostly this is done to allow 3 arguments in functions 
 *  that do not normally take them.
 *
 **************/

static R_INLINE double NORMAL           (double value, double *args, int nArgs){
    if(isnan(value)!=0) return 0;
    int useLowerTail = (round(pnorm(value, args[0], args[1], FALSE, TRUE))==0); 
    // Check for underflow.
    double funcVal= expDif(pnorm(value-0.5, args[0], args[1], 
        useLowerTail, TRUE), pnorm(value+0.5, args[0], args[1], 
        useLowerTail, TRUE));
    return funcVal; 
    }
static R_INLINE double dNORMAL      (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    return dnorm(value, args[0], args[1], TRUE);
    } // Sometimes we have continuous value ...
static R_INLINE double BETA         (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    int useLowerTail = (round(pbeta(value, args[0], args[1], FALSE, TRUE))==0); 
    // Check for underflow.
    return expDif(pbeta(value-0.5, args[0], args[1], useLowerTail, TRUE), 
        pbeta(value+0.5, args[0], args[1], useLowerTail, TRUE)); 
    }
static R_INLINE double NONCENTRALBETA   (double value, double *args, int nArgs){
    if(isnan(value)!=0) return 0;
    int useLowerTail = (round(pnbeta(value, args[0], args[1], args[2], 
        FALSE, TRUE))==0); // Check for underflow.
    return expDif(pnbeta(value-0.5, args[0], args[1], args[2], 
        useLowerTail, TRUE), pnbeta(value+0.5, args[0], args[1], args[2], 
        useLowerTail, TRUE)); 
    }
static R_INLINE double BINOMIAL     (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    return dbinom(value, args[0], args[1], TRUE);
    }
static R_INLINE double EXPONENTIAL  (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    int useLowerTail = (round(pexp(value, args[0], FALSE, TRUE))==0); 
    // Check for underflow.
    return expDif(pexp(value-0.5, args[0], useLowerTail, TRUE), 
        pexp(value+0.5, args[0], useLowerTail, TRUE)); 
    }
static R_INLINE double GAMMA          (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    int useLowerTail =(round(pgamma(value, args[0], args[1], FALSE, TRUE))==0);
    // Check for underflow.
    return expDif(pgamma(value-0.5, args[0], args[1], useLowerTail, TRUE), 
        pgamma(value+0.5, args[0], args[1], useLowerTail, TRUE)); 
    }
static R_INLINE double dGAMMA       (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    return dgamma(value, args[0], args[1], TRUE);
    } // Included for legacy support of the transcript detection HMM in 
    // Hah et al. 2011.
static R_INLINE double GAMMA_p1     (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    int useLowerTail = 
        (round(pgamma(value+1,args[0], args[1], FALSE, TRUE))==0);
    // Check for underflow.
    return expDif(pgamma(value-0.5+1, args[0], args[1], useLowerTail, TRUE), 
        pgamma(value+0.5+1, args[0], args[1], useLowerTail, TRUE)); 
    } // Gamma of random variable (x+1).
static R_INLINE double HYPERGEOMETRIC  (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    return dhyper(value, args[0], args[1], args[2], TRUE);
    }
static R_INLINE double NEGATIVEBINOMIAL (double value, double *args, int nArgs){
    if(isnan(value)!=0) return 0;
    return dnbinom(value, args[0], args[1], TRUE);
    }
static R_INLINE double POISSON      (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    return dpois(value, args[0], TRUE);
    }
static R_INLINE double UNIFORM      (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    int useLowerTail = (round(punif(value, args[0], args[1], FALSE, TRUE))==0); 
    // Check for underflow.
    return expDif(punif(value-0.5, args[0], args[1], useLowerTail, TRUE), 
        punif(value+0.5, args[0], args[1], useLowerTail, TRUE)); 
    }
// args order: alpha, mu, sigma, lambda
static R_INLINE double NORMAL_EXP     (double value, double *args, int nArgs) {
    if(isnan(value)!=0) return 0;
    // int useLowerTailN = (round(pnorm(value, args[1], 
    //  args[2], FALSE, TRUE))==0); // Check for underflow.  Unused
    // int useLowerTailE = (round(pexp(value, args[3], FALSE, TRUE))==0); 
    // Check for underflow. Unused
    double N= NORMAL(value, &args[1], 2);
    double E= EXPONENTIAL(value, &args[3], 1);
    double funcVal= expSum2(log(args[0])+N,log(1-args[0])+E);
    return funcVal;
    }
// Can also write a function for a multinomial distribution given an index of 
// possible emissions.

static R_INLINE double NORMAL_EXP_MINUS (double value, double *args, 
    int nArgs) {
    if(isnan(value)!=0) return 0;
    value= -1*value; // reverse the VALUE.
    // int useLowerTailN = 
    // (round(pnorm(value, args[1], args[2], FALSE, TRUE))==0); 
    // Check for underflow. Unused
    // int useLowerTailE = (round(pexp(value, args[3], FALSE, TRUE))==0); 
    // Check for underflow. Unused
    double N= NORMAL(value, &args[1], 2);
    double E= EXPONENTIAL(value, &args[3], 1);
    double funcVal= expSum2(log(args[0])+N,log(1-args[0])+E);
    return funcVal;
    }


#endif
