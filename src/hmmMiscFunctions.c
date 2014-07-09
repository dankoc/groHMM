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
 *  Misclaneous functions for the hmm implementation in the GROseq package ...
 *  Charles Danko.
 *
 *  Started 12-2-09
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
#include <assert.h>

#include "hmmHeader.h"

/*****************************************************************************
 *
 * Init functions for matrices.  R_alloc generally used ...
 *
 * matrix: rows x cols
 *   access m[row][col]
 *
 * if reverse = 1, cols are allocated from max(row) to 0, instead of from
 * 0 to max(row)
 *
 *****************************************************************************/
extern double **matrix_alloc(int rows, int cols, int reverse) {
  double **matrix = (double**)calloc(rows, sizeof(double*));
  int i;
  size_t block_len = cols * rows;
  /* allocate really large (hopefully contiguos) block */
  double *large_block = (double*)calloc(block_len, sizeof(double));
  double *col_ptr;

  if (reverse == 1) {
    for (i = rows - 1, col_ptr = large_block; i >= 0; --i, col_ptr += cols)
      matrix[i] = col_ptr;
  } else {
    for (i = 0, col_ptr = large_block; i < rows; ++i, col_ptr += cols)
      matrix[i] = col_ptr;
  }

  return matrix;
}

extern void matrix_free(double ** matrix, int reverse, int rows) {
  if (matrix == NULL)
    return;

  if (reverse == 1)
    free(matrix[rows - 1]);/* works because we allocated a single large block*/
  else
    free(matrix[0]); /* works because we allocated a single large block */
  free(matrix);
}


extern int **imatrix_alloc(int rows, int cols, int reverse) {
  int **matrix = (int**)calloc(rows, sizeof(int*));
  int i;
  size_t block_len = cols * rows;
  /* allocate really large (hopefully contiguos) block */
  int *large_block = (int*)calloc(block_len, sizeof(int));
  int *col_ptr;

  if (reverse == 1) {
    for (i = rows - 1, col_ptr = large_block; i >= 0; --i, col_ptr += cols)
      matrix[i] = col_ptr;
  } else {
    for (i = 0, col_ptr = large_block; i < rows; ++i, col_ptr += cols)
      matrix[i] = col_ptr;
  }

  return matrix;
}

extern void imatrix_free(int ** matrix, int reverse, int rows) {
  if (matrix == NULL)
    return;

  if (reverse == 1)
    free(matrix[rows - 1]);/* works because we allocated a single large block*/
  else
    free(matrix[0]); /* works because we allocated a single large block */
  free(matrix);
}



/**************
 * 
 *  Set up paremeters for an HMM.
 *
 *  Used to interface R and C in hmmViterbi.cpp and hmmEM.cpp.
 *
 **************/
extern hmm_t *setupHMM(SEXP nstates, SEXP emiprobDist, SEXP emiprobVars, 
    SEXP nEmis, SEXP tprob, SEXP iprob) {
    hmm_t *hmm = (hmm_t*)R_alloc(1, sizeof(hmm_t));

    hmm[0].n_states  = INTEGER(nstates)[0];
    hmm[0].n_emis    = INTEGER(nEmis)[0];
    hmm[0].log_iProb = REAL(iprob); // Initial probabilities.

    hmm[0].log_tProb = (double**)R_alloc(hmm[0].n_states, sizeof(double*)); 
        // (n_states x n_states).
    hmm[0].em_args   = (double**)
        R_alloc(hmm[0].n_states*hmm[0].n_emis, sizeof(double*)); 
        // (3 x n_states).

   /***************************************************
    * Set up the transition probability distributions.   
    ****************************************************/
    for(int s=0;s<(hmm[0].n_states);s++)
        hmm[0].log_tProb[s] = REAL(VECTOR_ELT(tprob, s));

   /***************************************************
    * Set up the emission probability distributions.   
    ****************************************************/
    for(int s=0;s<(hmm[0].n_states*hmm[0].n_emis);s++)
        hmm[0].em_args[s]   = REAL(VECTOR_ELT(emiprobVars, s));

    hmm[0].log_eProb = (emiss_func) R_alloc(hmm[0].n_states*hmm[0].n_emis, 
        sizeof( emiss_func ));
    for(int i=0;i<(hmm[0].n_states*hmm[0].n_emis);i++) {
        if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "norm") == 0)       
            hmm[0].log_eProb[i] = NORMAL;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "dnorm") == 0)     
            hmm[0].log_eProb[i] = dNORMAL;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "beta") == 0)      
            hmm[0].log_eProb[i] = BETA;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "nbeta") == 0)     
            hmm[0].log_eProb[i] = NONCENTRALBETA;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "binom") == 0)     
            hmm[0].log_eProb[i] = BINOMIAL;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "exp") == 0)       
            hmm[0].log_eProb[i] = EXPONENTIAL;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "gamma") == 0)     
            hmm[0].log_eProb[i] = GAMMA;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "dgamma") == 0)        
            hmm[0].log_eProb[i] = dGAMMA;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "gamma_scale1") == 0)  
            hmm[0].log_eProb[i] = GAMMA;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), 
            "gamma_SHAPEeq1overSCALE") == 0) 
            hmm[0].log_eProb[i] = GAMMA;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "gamma_p1") == 0)  
            hmm[0].log_eProb[i] = GAMMA_p1;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "hyper") == 0)     
            hmm[0].log_eProb[i] = HYPERGEOMETRIC;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "nbinom") == 0)    
            hmm[0].log_eProb[i] = NEGATIVEBINOMIAL;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "pois") == 0)      
            hmm[0].log_eProb[i] = POISSON;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "uniform") == 0)   
            hmm[0].log_eProb[i] = UNIFORM;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "normexp") == 0)       
            hmm[0].log_eProb[i] = NORMAL_EXP;
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "normexpminus") == 0)  
            hmm[0].log_eProb[i] = NORMAL_EXP_MINUS;
        else error("ERROR: Can't set up EM.  Emission distribution ('%s') \
            not recognized!", CHAR(STRING_ELT(emiprobDist, i)));
    }

    return(hmm);
}

/*******************************************************************
 *
 *  Update paremeters of the emission functions for EM ...
 *
 *  SHOULDN'T BE NECESSARY TO DO THIS ... POSTERIORS CAN BE USED DIRECTLY ...
 * 
 *******************************************************************/

/*
 *  Returns the sum of the elements of logValues.
 *  logValues -- double*, vector.
 *  length    -- number of elements in logValues.
 */

/////////////////////////////////////////////////////////////////////
// Gamma ////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

extern void* SSallocGamma(int num) {
  void* ss = (void*)calloc(1, sizeof(ssGamma)); 
  // Assumed to be ssGamma*; // Cast to void* ??
  ssGamma *SS = (ssGamma*)ss;
  SS[0].N          =0;
  SS[0].sumPiXi    =0;
  SS[0].sumPiXiSq  =0;
  SS[0].sumLogPiXi =0;
  return(ss);
}
extern void SStatsGamma(int state, int emis_indx, void* ss, fwbk_t fwbk) {
  double PP, logPP, data_i;
  double epsilon=0.001;
  ssGamma *SS = (ssGamma*)ss;
  for(int position=0;position<fwbk.N;position++) {
    if(isnan(fwbk.data[emis_indx][position]) != 0) continue;
    logPP =  MargainalizeSumLogProbOver(state, position, fwbk);

    // If the contribution matters at all within the bounds of the machine ... 
    // add its contribution.
    if(!(logPP <= epsilon)) Rprintf("[SSallocGamma] -- \
        Assertion about to fail!  logPP= %d\n", logPP);
    assert(logPP < epsilon); // Require to be less than rounding error...
    if(-1*logPP < APPROX_EXP_VALUE_THRESHOLD && 
        !isnan(fwbk.data[emis_indx][position])) {
        PP = exp(logPP);
        data_i = fwbk.data[emis_indx][position];
        if(data_i == 0) data_i = epsilon;
        SS[0].N          += PP;
        SS[0].sumPiXi    += PP*data_i;
        SS[0].sumPiXiSq  += PP*data_i*data_i;
        SS[0].sumLogPiXi += PP*log(data_i);
    }

  }
//  Rprintf("[SStatsGamma]: SS.N: %f; SS.sumPiXi: %f; SS.sumLogPiXi: %f\n", 
//  SS[0].N, SS[0].sumPiXi, SS[0].sumLogPiXi);
}
extern void SStatsGamma_p1(int state, int emis_indx, void* ss, fwbk_t fwbk) {
  double PP, logPP;
  ssGamma *SS = (ssGamma*)ss;
  for(int position=0;position<fwbk.N;position++) {
    logPP =  MargainalizeSumLogProbOver(state, position, fwbk);

    // If the contribution matters at all within the bounds of the machine ... 
    // add its contribution.
    assert(logPP <= 0);
    if(-1*logPP < APPROX_EXP_VALUE_THRESHOLD && 
        !isnan(fwbk.data[emis_indx][position])) {
        PP = exp(logPP);
        SS[0].N          += PP;
        SS[0].sumPiXi    += PP*(fwbk.data[emis_indx][position]+1);
        SS[0].sumLogPiXi += PP*log(fwbk.data[emis_indx][position]+1);
    }

  }
}
extern void UpdateGamma(int state, void* ss, hmm_t *hmm) {
  ssGamma *SS = (ssGamma*)ss;
  double *shape= (double*)Calloc(1, double);
  double *scale= (double*)Calloc(1, double);
  int updateRetVal= MLEGamma(SS[0].N, SS[0].sumPiXi, SS[0].sumLogPiXi, 
                        shape, scale);
  if(updateRetVal == 0) {
      hmm[0].em_args[state][0] = shape[0];
      hmm[0].em_args[state][1] = scale[0];
  }
  else {
    Rprintf("WARNING! [UpdateGamma]\t--> Gamma for state %d update failed \
        due to instibility!  Using Shape: %f; Scale: %f\n", 
        state, hmm[0].em_args[state][0], hmm[0].em_args[state][1]);
  }
  Free(shape); 
  Free(scale);
}
// Used to fit a constrained gamma, where E[x] = 1, and shape=1/scale.
extern void UpdateGamma_SHAPEeq1overSCALE(int state, void* ss, hmm_t *hmm) {
  ssGamma *SS = (ssGamma*)ss;
  MLEGamma_SHAPEeq1overSCALE(SS[0].N, SS[0].sumPiXi, SS[0].sumLogPiXi, 
    SS[0].sumPiXiSq, &(hmm[0].em_args[state][0]), &(hmm[0].em_args[state][1]));
  Rprintf("[UpdateGammaConstrained]\t--> Shape: %f; Scale: %f; \
    Shape/Scale: %f (shape/scale must be 1!)\n", hmm[0].em_args[state][0], 
    hmm[0].em_args[state][1], 
    (hmm[0].em_args[state][0]/hmm[0].em_args[state][1]));
}
// Used to fit a constrained gamma where E[mean] = E[var].
extern void UpdateGamma_SCALE1(int state, void* ss, hmm_t *hmm) {
  ssGamma *SS = (ssGamma*)ss;
  MLEGamma_SCALE1(SS[0].N, SS[0].sumPiXi, SS[0].sumLogPiXi, 
    &(hmm[0].em_args[state][0]), &(hmm[0].em_args[state][1]));
  Rprintf("[UpdateGamma_Scale1]\t--> Shape: %f; Scale: %f\n", 
    hmm[0].em_args[state][0], hmm[0].em_args[state][1]);
}
extern void SSfreeGamma(void* ss) { 
  free((ssGamma*)ss); 
}

/////////////////////////////////////////////////////////////////////
// Normal
/////////////////////////////////////////////////////////////////////

extern void* SSallocNormal(int num) {
  void* ss = (void*)calloc(1,sizeof(ssNormal));
  ssNormal *SS = (ssNormal*)ss;
  SS[0].N          = 0;
  SS[0].sumPiXi    = 0;
  SS[0].sumPiXiSq  = 0;
  return(ss);
}

extern void SStatsNormal(int state, int emis_indx, void* ss, fwbk_t fwbk) {
  double PP, logPP;
  ssNormal *SS = (ssNormal*)ss;
  for(int position=0;position<fwbk.N;position++) {
    if(isnan(fwbk.data[emis_indx][position]) != 0) continue;
    logPP =  MargainalizeSumLogProbOver(state, position, fwbk);

    // If the contribution matters at all within the bounds of the machine ... 
    // add its contribution.
    if(-1*logPP < APPROX_EXP_VALUE_THRESHOLD && 
        !isnan(fwbk.data[emis_indx][position])) {
        PP = exp(logPP);
        SS[0].N          += PP;
        SS[0].sumPiXi    += PP*fwbk.data[emis_indx][position];
        SS[0].sumPiXiSq  += 
            PP*fwbk.data[emis_indx][position]*fwbk.data[emis_indx][position];
    }
  }
}

extern void UpdateNormal(int state, void* ss, hmm_t *hmm) {
  ssNormal *SS = (ssNormal*)ss;
  double epsilon=0.00001;
  double *stateParams = hmm[0].em_args[state];

  // Update mean.  Mean is the first entry in the **em_args matrix.
  stateParams[0] = SS[0].sumPiXi/SS[0].N;

  // Update var= (Sum_i{Xi^2}/N) - (mean^2).
//  Rprintf("[UpdateNormal]\t--> N: %f; Mean: %f; sumPiXiSq: %f\n", 
//  SS[0].N, stateParams[0], SS[0].sumPiXiSq);
  stateParams[1] = SS[0].sumPiXiSq/SS[0].N-(stateParams[0]*stateParams[0]);

  assert(stateParams[1] > -1); // Rounding error can makes this quantity 
                // JUST BARELY <1 for samples with 0 varience.
                // This just checks that the error is not excessive 
                // (and thus explained by a bug in addition to rounding error).
  if(stateParams[1] < epsilon)
    stateParams[1] = epsilon; // Keep the varience >0.
  stateParams[1] = sqrt(stateParams[1]);

  Rprintf("[UpdateNormal]\t--> Mean: %f; Stdev: %f\n", 
    hmm[0].em_args[state][0], hmm[0].em_args[state][1]);
}
extern void SSfreeNormal(void* ss) { 
  free((ssNormal*)ss); 
}

/////////////////////////////////////////////////////////////////////
// Normal-Exp -- NOTE: DOSEN'T SUPPORT MISSING VALUES YET!
/////////////////////////////////////////////////////////////////////

extern void* SSallocNormExp(int num) {
  void* ss = (void*)calloc(1,sizeof(ssNormExp));
  ssNormExp *SS = (ssNormExp*)ss;
  SS[0].containsData = 0; // false;
  return(ss);
}
extern void SStatsNormExp(int state, int emis_indx, void* ss, fwbk_t fwbk) {
  double wi, *newEx;
  ssNormExp *SS = (ssNormExp*)ss;
  int arraySize=fwbk.N;
  int oldSize=0;

  // Count the number of positions and (re)-alloc ss.
  if(SS[0].containsData) {
    oldSize= (int)SS[0].ex[0];
    arraySize+= oldSize;
    newEx = (double*)calloc(arraySize*2+1, sizeof(double));
    for(int i_cp=1;i_cp<(2*oldSize+1);i_cp+=2)
      newEx[i_cp]=SS[0].ex[i_cp];
  }
  else { // Otherwise ... just alloc a new array.
    newEx = (double*)calloc(arraySize*2+1, sizeof(double));  
  }
  newEx[0] = (double)arraySize;
  
  // Put in new data from this chromsome.
  for(int position=0;position<fwbk.N;position++) {
    if(isnan(fwbk.data[emis_indx][position]) != 0) continue;
    newEx[(position+oldSize)*2+1] = fwbk.data[emis_indx][position]; // xi
    newEx[(position+oldSize)*2+2] = 
        exp(MargainalizeSumLogProbOver(state, position, fwbk)); // wi
  }
  if(SS[0].containsData)
     free((double*)SS[0].ex);
  SS[0].ex = newEx;
  SS[0].containsData = 1; //true;
}
extern void UpdateNormExp(int state, void* ss, hmm_t *hmm) {
  Rprintf("[UpdateNormExp] START");
  ssNormExp *SS = (ssNormExp*)ss;
  double epsilon=0.00001;
  double *stateParams = hmm[0].em_args[state];
  double *sp = (double*)calloc(4,sizeof(double));// Copy it...;
  double *fmin=(double*)calloc(1,sizeof(double));
  int *fail=(int*)calloc(1,sizeof(int));
  int MAXIT= 100;
  double tol=1e-3; // Don't have to be all that precise at each update ... 
                  // so long as we are making headway.
  double *ex= SS[0].ex; 
  int *fncount=(int*)calloc(1,sizeof(int)), 
    *grcount=(int*)calloc(1,sizeof(int));

  // Update mean.  Mean is the first entry in the **em_args matrix.
  cgmin(4, stateParams, sp, fmin, normal_exp_optimfn, 
    normal_exp_optimgr, fail, tol, tol, ex, 1, 0, fncount, grcount, MAXIT);

  if(fail[0] != 0) {
    Rprintf("[UpdateNormExp] WARNING::Updates failed w/ message %d.  \
        fncount= %d ; grcount= %d\n", fail[0], fncount[0], grcount[0]);
  }
  
  stateParams = sp;  
  // POSSIBLE MEMORY LEAK: Does this leak the old 'hmm[0].em_args[state]'??
  Rprintf("[UpdateNormExp]\t--> Alpha: %f; Mean: %f; Stdev: %f; Lambda: %f\n", 
    hmm[0].em_args[state][0], hmm[0].em_args[state][1], 
    hmm[0].em_args[state][2], hmm[0].em_args[state][3]);
}
extern void SSfreeNormExp(void* ss) { 
  ssNormExp *SS = (ssNormExp*)ss;
  free((double*)SS[0].ex);
  free((ssNormExp*)ss); 
}

/////////////////////////////////////////////////////////////////////
// Poisson
/////////////////////////////////////////////////////////////////////

extern void* SSallocPoisson(int num) { 
  void* ss = (void*)calloc(1,sizeof(ssPoisson)); 
  ssPoisson *SS = (ssPoisson*)ss;
  SS[0].N       =0;
  SS[0].sumPiXi =0;
  return(ss);
}
extern void SStatsPoisson(int state, int emis_indx, void* ss, fwbk_t fwbk) {
  double PP, logPP;
  ssPoisson *SS = (ssPoisson*)ss;
  for(int position=0;position<fwbk.N;position++) {
    if(isnan(fwbk.data[emis_indx][position]) != 0) continue;
    logPP =  MargainalizeSumLogProbOver(state, position, fwbk);

    if(-1*logPP < APPROX_EXP_VALUE_THRESHOLD && 
        !isnan(fwbk.data[emis_indx][position])) {
        PP = exp(logPP);
        SS[0].N          += PP;
        SS[0].sumPiXi    += PP*fwbk.data[emis_indx][position];
    }
  }

}
extern void UpdatePoisson(int state, void* ss, hmm_t *hmm) {
  ssPoisson *SS = (ssPoisson*)ss;
  hmm[0].em_args[state][0] = SS[0].sumPiXi/SS[0].N;
  Rprintf("[UpdatePoisson]\t--> Lambda: %f\n", hmm[0].em_args[state][0]);
}
extern void SSfreePoisson(void* ss) { 
  free((ssPoisson*)ss); 
} 

/////////////////////////////////////////////////////////////////////
// Exponential distribution 
/////////////////////////////////////////////////////////////////////

extern void *SSallocExp(int num) {
  void* ss = (void*)calloc(1,sizeof(ssExp)); 
  ssExp *SS = (ssExp*)ss;
  SS[0].N       =0;
  SS[0].sumPiXi =0;
  return(ss);
}
extern void SStatsExp(int state, int emis_indx, void* ss, fwbk_t fwbk) {
  double PP, logPP;
  ssExp *SS = (ssExp*)ss;
  for(int position=0;position<fwbk.N;position++) {
    if(isnan(fwbk.data[emis_indx][position]) != 0) continue;
    logPP =  MargainalizeSumLogProbOver(state, position, fwbk);

    if(-1*logPP < APPROX_EXP_VALUE_THRESHOLD && 
        !isnan(fwbk.data[emis_indx][position])) {
        PP = exp(logPP);
        SS[0].N          += PP;
        SS[0].sumPiXi    += PP*fwbk.data[emis_indx][position];
    }
  }
}
extern void UpdateExp(int state, void* ss, hmm_t *hmm) {
  ssExp *SS = (ssExp*)ss;
  hmm[0].em_args[state][0] = SS[0].N/SS[0].sumPiXi; // Lambda = 1/mean.
  Rprintf("[UpdateExp]\t--> Lambda: %f\n", hmm[0].em_args[state][0]);
}
extern void SSfreeExp(void* ss) {
  free((ssExp*)ss); 
}

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
// Transition probibililty update -- from state k to state l 
// num --> the number of states in the model.
// sequences --> the number of training sequences.
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

extern void* TransAlloc(int num, int sequences) { 
  void* ss = (void*)calloc(1,sizeof(ssTransition)); 
  ssTransition *SS = (ssTransition*)ss;
  SS[0].totalTransK = matrix_alloc(num, sequences, FALSE);
  return(ss);
}

extern void  TransUpdate(int state, int sequence, void* ss, emiss_func EMI, 
    fwbk_t fwbk) {
  ssTransition *SS = (ssTransition*)ss;
  int N = fwbk.N;
  int n_states = fwbk.hmm[0].n_states;
  int n_emis = fwbk.hmm[0].n_emis;
  int emis_count;
  double CurrentSum, ChromSum;
  double scalefactor;

  for(int l=0;l<n_states;l++) {
    double *A = (double*)calloc(N, sizeof(double));
    ChromSum =0;

    /* Finds a suitable maximum for this chromosome ... */
    A[0] = fwbk.forward[0][state]+fwbk.backward[1][l]+
        fwbk.hmm[0].log_tProb[state][l];//+
    for(emis_count=0;emis_count<n_emis;emis_count++) 
        A[0] += (EMI[l+n_states*emis_count])(fwbk.data[emis_count][1], 
            fwbk.hmm[0].em_args[l+n_states*emis_count], 4);

    scalefactor= A[0];
    for(int i=1;i<(N-1);i++) {
        A[i] = fwbk.forward[i][state]+fwbk.backward[i+1][l]+
            fwbk.hmm[0].log_tProb[state][l];//+
        for(emis_count=0;emis_count<n_emis;emis_count++)
            A[i] += (EMI[l+n_states*emis_count])(fwbk.data[emis_count][i+1], 
                fwbk.hmm[0].em_args[l+n_states*emis_count], 4);

        scalefactor = max(scalefactor, A[i]);
    }

    /* Calculates sum over chromsome using scalefactor ... */
    for(int i=0;i<(N-1);i++) {
        CurrentSum = A[i]-scalefactor;
        if(-1*CurrentSum < APPROX_EXP_VALUE_THRESHOLD)
            ChromSum += exp(CurrentSum);
    }
    free(A);

    assert(!isnan(ChromSum)); // Try this instead ...
    
    /* Assign the sum to the proper chromosome in total transitions. */
    SS[0].totalTransK[l][sequence] = log(ChromSum) + scalefactor - fwbk.log_px;

    Rprintf("[TransUpdate]\t--> Chrom: %d; State: %d; ChromSum=%f; Final=%f\n", 
            sequence, l, ChromSum, SS[0].totalTransK[l][sequence]);
  }
}

extern void  TransUpdateP(int state, int nSequences, void* ss, hmm_t *hmm) {
  ssTransition *SS = (ssTransition*)ss;

  // Calculate Akl ==> i.e. the total times in state k; sum over other states, 
  // l' of: {k->l' transition}.
  // Durbin 3.20 sum over j's (training sequences) for each k->l transition.
  double *ExpectedTransitions = 
    (double*)calloc(hmm[0].n_states, sizeof(double));
  for(int l=0;l<hmm[0].n_states;l++){
    ExpectedTransitions[l] = expSum(SS[0].totalTransK[l], nSequences); 
    // Foreach state, state->l, sum over training sequences.
  }

  // Calculate sum_l' {Akl'}  --> i.e. Bottom of equation 3.18.
  double TotalSum = expSum(ExpectedTransitions, hmm[0].n_states);

  // Update paremeters for each state ...
  for(int l=0;l<hmm[0].n_states;l++) {
    double CurrentDiff = ExpectedTransitions[l] - TotalSum;

    // Error checking.
    if(isnan(CurrentDiff)) {
        Rprintf("ASSERTION ABOUT TO FAIL.  CurrentDiff= %f\n", CurrentDiff);
        error("CurrentDiff is nan.");
    }
    assert(CurrentDiff <= 0.001); // Allow for some rounding error...

    hmm[0].log_tProb[state][l] = CurrentDiff;

    Rprintf("[UpdateTransitionProb]\t--> TP_{%d->%d}: %f\n", state, 
        l, hmm[0].log_tProb[state][l]);
  }

  free(ExpectedTransitions);
}

extern void  TransFree(void* ss) {
  ssTransition *SS = (ssTransition*)ss;
  matrix_free(SS[0].totalTransK, FALSE, -1);
  free((ssTransition*)ss);
}

