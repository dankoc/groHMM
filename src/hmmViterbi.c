/***************************************************************************
**
**   Copyright 2009, 2010, 2011 Charles Danko and Andre Martins.
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
 *  Source code written for GRO-seq package by Charles Danko.
 *
 *  Implementation of standard Viterbi algorithm heavily based off of 
 *  code by Andre Martins (alm253@cornell.edu).
 *
 *  This implementation assumes that emission probabilities are modeled by 
 *  parametric probability distributions.  
 *  
 *  It should be fairly easy to write your own distribution/ function.  
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

/*******************************************************************************
 *
 * Viterbi Algorithm.  Initial implementation by Andre Martins.  
 * Motified by Charles.
 *
 * hmm - hmm structure with appropriate functions
 * data - user data pointer
 * seq_len - sequence length
 * matrix - (seq_len x hmm->n_states) preallocated matrix
 *          if NULL, viterbi will allocate internally
 * backptrs - (seq_len x hmm->n_states) preallocated matrix
 *            if NULL, viterbi will allocate internally
 *
 * Returns vector of size seq_len with the state path (excluding start state)
 *
 ******************************************************************************/
extern void viterbi_path(hmm_t hmm, double **data, int seq_len, 
    double **matrix, int **backptr, int *path) {
  int cols = hmm.n_states;
  int rows = seq_len;
  int m_free = (matrix == NULL);
  int b_free = (backptr == NULL);
  double  *log_iProb   = hmm.log_iProb;
  double **log_tProb   = hmm.log_tProb;
  emiss_func log_eProb = hmm.log_eProb;
  double **emargs      = hmm.em_args;
  int n_emis           = hmm.n_emis;

  int i, k, l, emis_count;
  double * m_row, * m_row_prev;
  int * b_row;
  int z;

  /* setup matrices */
  if (matrix == NULL)
    matrix = matrix_alloc(rows, cols, 0);

  if (backptr == NULL) {
    backptr = imatrix_alloc(rows, cols, 0);
  }

  /* fill first row */
  for (l = 0; l < cols; ++l) { /* -1 = start state */
    matrix[0][l] = log_iProb[l];
    for(emis_count=0;emis_count<n_emis;emis_count++) 
        matrix[0][l] += (log_eProb[l+cols*emis_count])(data[emis_count][0], 
            emargs[l+cols*emis_count], 4);
    backptr[0][l] = -1; /* stop */
  }

  /* inner rows */
  for (i = 1; i < rows; ++i) {
    m_row = matrix[i];
    m_row_prev = matrix[i - 1];
    b_row = backptr[i];

    for (l = 0; l < cols; ++l) {
        double max = -HUGE_VAL;
        int argmax = -1;
        
        for (k = 0; k<hmm.n_states; k++) {
            double value = m_row_prev[k] + log_tProb[k][l];

            if (value > max) {
                max = value;
                argmax = k;
            }
        }

        // for loop over emissions...
        m_row[l] = max; 
        for(emis_count=0;emis_count<n_emis;emis_count++) 
                m_row[l] += (log_eProb[l+cols*emis_count])(data[emis_count][i], 
                    emargs[l+cols*emis_count], 4);
        b_row[l] = argmax;
    }
  }
  
  /* backtrace */
  /* last state */
  i = rows - 1;
  double max = -HUGE_VAL;
  int argmax = -1;
  m_row = matrix[i];

  for (k = 0; k < cols; ++k) {
    double value = m_row[k];
    if (value > max) {
      max = value;
      argmax = k;
    }
  }
  path[i] = argmax;
  
  /* other states */
  z = path[i];
  for (l = rows - 1, i = rows - 2; l > 0; --l, --i) {
    z = backptr[l][z];
    path[i] = z;
  }

  // Free matrices.
  if(m_free)
      matrix_free(matrix, 0, -1);
  if(b_free)
    imatrix_free(backptr, 0, -1);
}

/*******************************************************************************
 *
 * Rviterbi -- Returns a vector of states which maxamizes the probability given 
 *  the emissions. Serves as a wrapper to Andre's Viterbi implementation.
 *
 *  emi --  Vector of observed emission over all sequence.
 *  nEmis   --  Number of emissions vectors.
 *  nstates --  Number of states in the HMM.
 *  emiprobD--  Vector (1 x nstates) of strings representing the emission 
 *      probability distirubitons.
 *  emiprobV--  Matrix (3? x nstates) representing arguments to pdist*.  
 *      Unused paremeters can be set to 0.
 *  tprob   --  Matrix (nstates x nstates) representing transition 
 *      probabilities between states.
 *  iprob   --  Vector (1 x nstates) of initial probabilities.
 *
 *  Assumptions: 
 *  (1) Emission probabilities are drawn from a continuous prob. distribution;
 *      which prob. distribution they are drawn from is set using emiprobD.
 *  (2) emi represents emissions from one sequence only!  Rviterbi, and/or 
 *      viterbi_path are run multiple times, once for each sequence.
 *
 *  Questions/Implementation Notes:
 *  What is the best way to pass a distribution from R?  
 *    Can I pass the name of the function?!
 *    Also needs to pass the arguments!
 *    For now, assume that we can pass a string representing the function name, 
 *    and make it in an if statement.
 *
 *  Updates/Information:
 *  2009-09-09: Wrote this function to integrate R with Andre's 
 *  Viterbi implementation.
 *
 ******************************************************************************/
SEXP Rviterbi(SEXP emi, SEXP nEmis, SEXP nstates, SEXP emiprobDist, 
    SEXP emiprobVars, SEXP tprob, SEXP iprob) {

    // Init hmm_t.
    hmm_t *hmm = setupHMM(nstates, emiprobDist, emiprobVars, nEmis, tprob, 
        iprob);

    // Set up emissions.
    int maxT = Rf_nrows(VECTOR_ELT(emi, 0));
    double **emisDATA = (double**)R_alloc(hmm[0].n_emis, sizeof(double*));
    for(int s=0;s<hmm[0].n_emis;s++) {
        emisDATA[s]   = REAL(VECTOR_ELT(emi, s));
    }
    
    SEXP hiddenStatesR;
    PROTECT(hiddenStatesR = allocVector(INTSXP,maxT));
    int *hiddenstates = INTEGER(hiddenStatesR);
    viterbi_path(hmm[0], emisDATA, maxT, NULL, NULL, hiddenstates);

    // Get it in the correct state to return to R.
    UNPROTECT(1);

    // Return to R and exit!
    return(hiddenStatesR);
}

