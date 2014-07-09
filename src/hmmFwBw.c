/***************************************************************************
**
**   Copyright 2009, 2010, 2011 Andre Martins and Charles Danko.
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
 *  Forward/Backward algorithms.  Initial implementation by Andre Martins.
 *
 *  2009-11-23 File pulled in from implementation by Andre Martins (thanks!). 
 *
 *  TODO: 
 *  (1) Finish up backwards.
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

/*************************************************************************
 *
 * Convenient forward/backward data storage struct, alloc, and free functions.
 *
 * N        - sequence size
 * data     - observed emissions 
 * forward  - N X hmmt->n_states matrix representing the forward algorithm 
 *  output
 * backward - ^-- same; for backward algorithm
 * log_px   - log probability of model|data (???)
 * bk_log_px    - (???)
 *
 *************************************************************************/
extern fwbk_t * fwbk_alloc(double **data, int N, hmm_t *hmm) {
  fwbk_t * res = (fwbk_t*) calloc(1, sizeof(fwbk_t));
  res->forward = matrix_alloc(N, hmm->n_states, 0);
  res->backward = matrix_alloc(N, hmm->n_states, 1); 
  // Allocate this backwards, but then use it backwards?
  res->N = N;
  res->log_px = 0.0;
  res->bk_log_px = 0.0;
  res->data = data;
  res->hmm = hmm;
  return res;
}

extern void fwbk_free(fwbk_t * data) {
  matrix_free(data->forward, 0, data->N);
  matrix_free(data->backward, 1, data->N);
  free(data);
}


/******************************************************************
 *
 *  Thanks to Andre Martins for these forward/backward implementaitons!
 *
 *  Forward/Backward algorithms assumptions:
 *  - transition/emission functions return log's of probabilities
 *  - number of states does not include initial state
 *  - there is no end state
 *  - states numbered from 0 to (n_states-1)
 *  - data has been pre-allocated (N by n_states matrices)
 *
 ******************************************************************/
void forward(fwbk_t *data) {
// Vars for the implementation.
  int i, l, k, emis_count; // i --> position in sequence; l,k --> current state.
  double *m_col, *m_col_prev;
  double sum; // Note that sum is NOT in log-space.

// Shorthand for vars stored in hmm_t and fwbk_t structs.
  int n                = data->hmm->n_states;
  int n_emis           = data->hmm->n_emis;
  double  *log_iProb   = data->hmm->log_iProb;
  double **log_tProb   = data->hmm->log_tProb;
  emiss_func log_eProb = data->hmm->log_eProb;
  double **emargs    = data->hmm->em_args;

  int N                = data->N;
  double **matrix      = data->forward;
  double  **edata       = data->data; /* external data */
  double scalefactor, current_sum;

// Initialization step @ i=0.
  m_col = matrix[0];
  for (k = 0; k < n; ++k) {
    m_col[k] =  log_iProb[k]; 
    //    + (log_eProb[k])(edata[0], emisargs[k][0], emisargs[k][1], 
    //    emisargs[k][2]);
    for(emis_count=0;emis_count<n_emis;emis_count++) 
        matrix[0][k] += (log_eProb[k+n*emis_count])(edata[emis_count][0], 
            emargs[k+n*emis_count], 4);
  }
// Recursion step.
  for (i = 1; i < N; ++i) {
    m_col_prev = matrix[i - 1];
    m_col = matrix[i];

    for (l = 0; l < n; ++l) { // foreach state.

      // Compute inner sum across states.
      sum = 0;

      // Modified 12-4 to prevent underflow.
      // Calculating exp<-735.-somthing gives underflow.  
      // The solution works because: log(exp(a)+exp(b)) = 
      //    log(exp(a+n)+exp(b+n))-n
      // Choose n=-max(val in sum); I should likely check to make sure that 
      //    this is defined?!
      // Foreach value in the sum, check that n-valInSum does not 
      //    over/under-flow.
      // n ==> scalefactor.
      // Thanks, Melissa Hubisz!  Also mentioned in Durbin et. al.'s HMM book.
      scalefactor = (m_col_prev[0] + log_tProb[0][l]); // Init to first value.
      for(k=1; k<n; k++)
        scalefactor = max((m_col_prev[k] + log_tProb[k][l]), scalefactor);

        for (k = 0; k<n; k++) {
        current_sum = m_col_prev[k] + log_tProb[k][l] - scalefactor;

//  The following assertion fails if k --> l is prohibited by transition 
//      probabilities.
//       Happens because (scalefactor --> -inf) (it is the max); 
//   This causes (currentsum --> nan) because -inf+inf --> nan
//  if(!(current_sum <= 0))
//        assert(current_sum <= 0);
        if((-1*current_sum) < APPROX_EXP_VALUE_THRESHOLD)
          sum += exp(current_sum);

        if(i>(N-2) || i<2) // Report the first and last case for debuging...
            Rprintf("i=%d, l=%d, k=%d, prev[k]=%f, scalefactor=%f, \
                prod=%f, sum=%f\n", i, l, k, m_col_prev[k], 
                scalefactor, current_sum, sum);
      }

      // Update matrix.
      m_col[l] = log(sum)+scalefactor;
      for(emis_count=0;emis_count<n_emis;emis_count++) {
        m_col[l] += (log_eProb[l+n*emis_count])(edata[emis_count][i], 
            emargs[l+n*emis_count], 4);
      }
    }
  }

// Termination; assume no last state, so no transition prob ?! */
  m_col = matrix[N-1];
  sum=0;

  scalefactor = (m_col[0]);
  for (i = 1; i < n; ++i)
    scalefactor = max(m_col[i], scalefactor);

  for (i = 0; i < n; ++i) {
    current_sum = m_col[i]-scalefactor;
    if(!(current_sum <= 0)) {
      Rprintf("WARNING: Assertion about to fail in hmmFwBw.cpp (at line ~189).\
        current_sum= %f, m_col[%d]= %f, scalefactor= %f\n", 
        current_sum, i, m_col[i], scalefactor); //likely nan
      error("ERROR: current_sum <= 0 (likely NaN)\n"); 
      // If this fails, likely one of yoru stats is unreachable. 
      // Commented Apr. 10 2010.  Getting ALL rate genes (for NH and LC reps) 
      // to run through.
    }

    if(-1*(current_sum) < APPROX_EXP_VALUE_THRESHOLD)
      sum += exp(current_sum);
  }

  data->log_px = log(sum)+scalefactor;
//  Rprintf("[forward] Probability: %f\n", data->log_px);
}

/********************************
 * Backward algorithm.
 ********************************/
void backward(fwbk_t *data) {
  int i, l, k, emis_count, sf;
  double * m_col, * m_col_next;
  double sum; // Again, not in log-space.

// Shorthand for vars stored in hmm_t and fwbk_t structs.
  int n                = data->hmm->n_states;
  int n_emis           = data->hmm->n_emis;
  double  *log_iProb   = data->hmm->log_iProb;
  double **log_tProb   = data->hmm->log_tProb;
  emiss_func log_eProb = data->hmm->log_eProb;
  double **emargs    = data->hmm->em_args;

  int N                = data->N;
  double **edata        = data->data; // external data.
  double ** matrix     = data->backward;
  double scalefactor, current_sum;

  // Border conditions.
  m_col = matrix[N-1];
  for (k = 0; k < n; ++k)
    m_col[k] = 0; /* log(1) */

  // Inner cells.
  for (i = N - 2; i >= 0; --i) {
    m_col = matrix[i];
    m_col_next = matrix[i + 1];

    for (k = 0; k < n; ++k) {
      // Compute inner sum.
      sum=0;
      scalefactor = m_col_next[0] + log_tProb[k][0];
      for(emis_count=0;emis_count<n_emis;emis_count++) 
        scalefactor += (log_eProb[0+n*emis_count])
            (edata[emis_count][i+1], emargs[0+n*emis_count], 4);

      for (l=1;l<n;l++) {
        sf =  m_col_next[l] + log_tProb[k][l];
        for(emis_count=0;emis_count<n_emis;emis_count++) 
            sf += (log_eProb[l+n*emis_count])
                (edata[emis_count][i+1], emargs[l+n*emis_count], 4);
        scalefactor = max(scalefactor, sf);
      }

      for (l=0;l<n;l++) {
        current_sum = m_col_next[l] + log_tProb[k][l] -scalefactor;
        // +  // Are these indices correct ??? (99% certain)
        for(emis_count=0;emis_count<n_emis;emis_count++) 
            current_sum += (log_eProb[l+n*emis_count])(edata[emis_count][i+1], 
                emargs[l+n*emis_count], 4);
        if((-1*current_sum) < APPROX_EXP_VALUE_THRESHOLD)
          sum += exp(current_sum);
      }

      m_col[k] = log(sum) + scalefactor;
    }
  }

  // Termination. Calculate backward log-prob.
  m_col = matrix[0];
  sum=0;
  scalefactor = (m_col[0] + log_iProb[0]);
    // + (log_eProb[0])(edata[0], emisargs[0][0], emisargs[0][1], 
    // emisargs[0][2]));
  for(emis_count=0;emis_count<n_emis;emis_count++) 
    scalefactor += (log_eProb[0+n*emis_count])(edata[emis_count][0], 
        emargs[0+n*emis_count], 4);

  for (k = 1; k < n; ++k) { // BUG?!  Should be log_iProb[k]?! 
    sf =  m_col[k] + log_iProb[k];
    for(emis_count=0;emis_count<n_emis;emis_count++) 
        sf += (log_eProb[k+n*emis_count])(edata[emis_count][0], 
            emargs[k+n*emis_count], 4);
    scalefactor = max(scalefactor, sf);
  }

  for(k=0;k<n;k++) {
    current_sum = m_col[k]+ log_iProb[k]- scalefactor;
    for(emis_count=0;emis_count<n_emis;emis_count++) 
        current_sum += (log_eProb[k+n*emis_count])(edata[emis_count][0], 
        emargs[k+n*emis_count], 4);

    if((-1*current_sum) < APPROX_EXP_VALUE_THRESHOLD)
      sum += exp(current_sum);// Are these indices correct ???
  }
  data->bk_log_px = log(sum)+scalefactor;
}

