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
 *  Baum Welch EM implementation -- Written for the GRO-seq package by 
 *  Charles Danko.
 *
 *  TODO: 
 *    (1) Look into R support for multi-threading.
 *    (2) ...
 *
 *  2009-11-23 Imported and revised Forward/backward implementation 
 *      by Andre Martins.
 *  2009-11-24 Implemented Baum Welch algorthm.
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


/* Sets up the sufficient stats void* data types for use during EM. */
em_t *setupEM(hmm_t *hmm, SEXP emiprobDist, SEXP updatetrans, SEXP updateemis) {
    em_t *em = (em_t*)R_alloc(1, sizeof(em_t));

    /***************************************************
    * Set up the transition probability distributions.   
    ****************************************************/
    em->AllocTssFunc = (alloc_trans_sstats) 
        R_alloc(hmm->n_states, sizeof( alloc_trans_sstats ));
    em->UpdateTssFunc = (update_trans_SS) 
        R_alloc(hmm->n_states, sizeof( update_trans_SS ));
    em->UpdateTPFunc = (update_trans_Prob)  
        R_alloc(hmm->n_states, sizeof( update_trans_Prob ));
    em->FreeTssFunc = (free_trans_sstats)  
        R_alloc(hmm->n_states, sizeof( free_trans_sstats ));
    em->TransSS = (void**)R_alloc(hmm->n_states, sizeof(void*));

    for(int i=0;i<(hmm->n_states);i++) {
      em->AllocTssFunc[i]   = TransAlloc;
      em->UpdateTssFunc[i]  = TransUpdate;
      em->UpdateTPFunc[i]   = TransUpdateP;
      em->FreeTssFunc[i]    = TransFree;
    }

    /***************************************************
    * Set up the emission probability distributions.   
    ****************************************************/
    em->sstats_alloc = (alloc_emis_sstats) 
        R_alloc(hmm->n_states*hmm->n_emis, sizeof( alloc_emis_sstats ));
    em->sstats_emis  = (update_sstat_func) 
        R_alloc(hmm->n_states*hmm->n_emis, sizeof( update_sstat_func ));
    em->update_emis  = (update_emiss_func) 
        R_alloc(hmm->n_states*hmm->n_emis, sizeof( update_emiss_func ));
    em->free_emis_s  =  (free_emis_sstats) 
        R_alloc(hmm->n_states*hmm->n_emis, sizeof( free_emis_sstats  ));
    em->ss = (void**)R_alloc(hmm->n_states, sizeof(void*));

    for(int i=0;i<(hmm->n_states*hmm->n_emis);i++) {
        if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "norm") == 0 ||
           strcmp(CHAR(STRING_ELT(emiprobDist, i)), "dnorm") == 0 ) {
            em->sstats_alloc[i] = SSallocNormal;
            em->sstats_emis[i]  = SStatsNormal;
            em->update_emis[i]  = UpdateNormal;
            em->free_emis_s[i]  = SSfreeNormal;
        }
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "gamma") == 0 ||
                strcmp(CHAR(STRING_ELT(emiprobDist, i)), "dgamma") == 0 )   {
            em->sstats_alloc[i] = SSallocGamma;
            em->sstats_emis[i]  = SStatsGamma;
            em->update_emis[i]  = UpdateGamma;
            em->free_emis_s[i]  = SSfreeGamma;
        }
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "normexp") == 0 || 
                strcmp(CHAR(STRING_ELT(emiprobDist, i)), "normexpminus") == 0){
            em->sstats_alloc[i] = SSallocNormExp;
            em->sstats_emis[i]  = SStatsNormExp;
            em->update_emis[i]  = UpdateNormExp;
            em->free_emis_s[i]  = SSfreeNormExp;
        }
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "pois") == 0)  {
            em->sstats_alloc[i] = SSallocPoisson;
            em->sstats_emis[i]  = SStatsPoisson;
            em->update_emis[i]  = UpdatePoisson;
            em->free_emis_s[i]  = SSfreePoisson;        
        }
        else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "exp") == 0)   {
            em->sstats_alloc[i] = SSallocExp;
            em->sstats_emis[i]  = SStatsExp;
            em->update_emis[i]  = UpdateExp;
            em->free_emis_s[i]  = SSfreeExp;        
        }
        else error("Distribution ('%s') not recognized!", 
            CHAR(STRING_ELT(emiprobDist, i)));
    }
    
    /* booleans indicating whether or not to update by EM */
    em->updateTrans = INTEGER(updatetrans);
    em->updateEmis  = INTEGER(updateemis);

    return(em);
}
/* Parses C data types back into a format that can be easily imported 
 * back into the R session */
SEXP getEMReturnRTypes(hmm_t *hmm, int n_seq, SEXP emiprobVars, SEXP tprob, 
    SEXP emi, SEXP output) {
    int outopt = INTEGER(output)[0];
    fwbk_t *fwbk;
    double Q;

    // Add final model paremeters to the list.
    int RETURN_SIZE = 3, RETURN_INDX = 0;
    if(outopt>1) RETURN_SIZE++;
    if(outopt==10) RETURN_SIZE++;

    SEXP ListObject, hiddenStatesR, posteriors, postTrans; 
    // ListObject can be named more easily in the R function that calls it.
    PROTECT(ListObject = allocVector(VECSXP, RETURN_SIZE));
    SET_VECTOR_ELT(ListObject, RETURN_INDX++, emiprobVars);
    SET_VECTOR_ELT(ListObject, RETURN_INDX++, tprob);

    // Run Viterbi to get the most likely state paths?!
    SET_VECTOR_ELT(ListObject, RETURN_INDX++, 
        hiddenStatesR=allocVector(VECSXP, n_seq));
    if(outopt > 1) {
        SET_VECTOR_ELT(ListObject, RETURN_INDX++, 
            posteriors=allocVector(VECSXP, n_seq));
    }
    if(outopt == 10) {
        SET_VECTOR_ELT(ListObject, RETURN_INDX++, 
            postTrans=allocVector(VECSXP, n_seq));
    }

    for(int seq=0;seq<n_seq;seq++) {

        // Data in types.
        int maxT = Rf_nrows(VECTOR_ELT(emi, seq));
        double **data = (double**)R_alloc(hmm->n_emis, sizeof(double*));
        for(int s=0;s<hmm->n_emis;s++) {
            data[s]   = REAL(VECTOR_ELT(emi, seq+s*n_seq)); 
            //Correctly indexed?!
            assert(Rf_nrows(VECTOR_ELT(emi, seq+s*n_seq))==maxT); // Make sure!
        }
        
        // Return types.
        SET_VECTOR_ELT(hiddenStatesR, seq, allocVector(INTSXP,maxT));
        int *hiddenstates = INTEGER(VECTOR_ELT(hiddenStatesR, seq));

        viterbi_path(hmm[0], data, maxT, NULL, NULL, hiddenstates);

        // CGD --> 6-6-2011 --> Optionally also returns posteriors 
        // (only one sequence currently supported).
        if(outopt > 1) {
            // Set up R data types.  
            // How to set up a matrix?!
            SET_VECTOR_ELT(posteriors, seq, 
                allocMatrix(REALSXP, hmm->n_states, maxT)); 
                // REALSXP an actual type?!
            double *Rpost = REAL(VECTOR_ELT(posteriors, seq));
            
            // Run forward/ backward.
            fwbk= fwbk_alloc(data /* Emission data */, 
                    maxT /*(size of *data)*/, hmm);
            forward(fwbk);  // (2a)
            backward(fwbk); // (2a)
            Q= fwbk[0].log_px;

            // Foreach state, calculate posteriors.  Construction similar to: 
            // http://stat.ethz.ch/R-manual/R-devel/doc/manual/
            // R-exts.html#Attributes
            for(int state=0;state<hmm->n_states;state++)
              for(int position=0;position<maxT;position++)
                Rpost[state + hmm->n_states*position] = 
                    fwbk[0].forward[position][state]+
                        fwbk[0].backward[position][state]-Q; 
                        // Durbin book, eq. 3.14.
                
            if(outopt == 10) { 
                // Return post. prob. of transition between states 2 and 3.
                SET_VECTOR_ELT(postTrans, seq, allocVector(REALSXP, (maxT-1)));
                double *post_trans = REAL(VECTOR_ELT(postTrans, seq));
                for(int position=0;position<(maxT-1);position++) {
                  post_trans[position] = 
                    fwbk[0].forward[position][1]
                        +fwbk[0].backward[position+1][2]
                        +hmm->log_tProb[1][2] -Q; 
                        // +emission // Durbin book, eq. 3.19.
                  for(int emis_count=0;emis_count<hmm->n_emis;emis_count++) 
                    post_trans[position] += 
                        (hmm->log_eProb[2+hmm->n_emis*emis_count])
                        (data[emis_count][position+1], 
                        hmm->em_args[2+hmm->n_emis*emis_count], 4);
                }
            }
            
            // Cleanup...
            fwbk_free(fwbk);
       }
    }
    
    // Return list object.
    UNPROTECT(1);
    return(ListObject);
}

/*******************************************************************************
 *
 * rBaumWelchEM -- Implementation of the Baum Welch Expectation Maximazation 
 * algorithm used to estimate the paremeters of an HMM.
 *
 *  nstates     --  Number of states in the HMM.
 *  emi         --  List of real vectors representing the emissions over 
 *      different 
 *      sequences (e.g. chromosomes).
 *  nEmis       --  Number of emission probabilities for each sequence.
 *  emiprobD    --  Vector (1 x nstates) of strings representing the 
 *      emission probability distirubitons.
 *  emiprobV    --  Matrix (3 x nstates) representing arguments to pdist*.  
 *      Unused paremeters can be set to 0.
 *  tprob       --  Matrix (nstates x nstates) representing transition 
 *      probabilities between states.
 *  iprob       --  Vector (1 x nstates) of initial probabilities.
 *  threshold   --  Real. Threshold for convergence.
 *  updateTrans --  Boolean array.  Rows in the transition parameter 
 *      matrix that should be updated by EM.
 *  updateEmis  --  Boolean... [TODO: This should be bool* -- which state's 
 *      params are updated.]
 *  output      --  Integer.  Controls the output list options.
 *              =0|1 -- Standard, data+viterbi.
 *              >1 -- Returns posteriors.
 *              =10-- Returns posteriors & posterior prob. of transition 
 *                  from state 2->3.
 *  verbose     --  Integer. Controls the amount of output text.  
 *              =0 -- quiet, no output; 
 *              =1 -- regular updates;
 *
 *  Questions/Implementation Notes:
 *  -- Can be used to update transition probabilities, emission probabilities, 
 *      or both using the update variables.
 *  -- Currently no support for using multiple threads on multithread 
 *      processors.  This may come in the future depending on how support 
 *      is done in R.  It is also possible to implement this support
 *      on the C-side, but will require additional dependencies.
 *  -- Removes data from each hmm run as quickly as possible to prevent 
 *      memory from blowing up!
 *  -- Nasty long function.  Considering re-factoring R-side into separate 
 *      functions (e.g. R_var_setup(); R_return_var()).
 *  
 *  Updates/Information:
 *  2009-11-24 Wrote this function.
 *
 ******************************************************************************/
SEXP RBaumWelchEM(SEXP nstates, SEXP emi, SEXP nEmis, SEXP emiprobDist, 
    SEXP emiprobVars, SEXP tprob, SEXP iprob, SEXP threshold, SEXP updatetrans, 
    SEXP updateemis, SEXP output, SEXP verbose) {

/*************************************************************
 *  Set up data structures that are used in EM.
 *************************************************************/

    /* Init logical variables. */
    int verb = INTEGER(verbose)[0];
    
    if(verb) Rprintf("Initializing Baum-Welch EM.\n");

    /* likelihoods */
    double Qprev = -999999999999999.0;
    double Q, Qdiff;
    double T = REAL(threshold)[0];

    /* transfer info from R vars into hmm struct */
    hmm_t *hmm = setupHMM(nstates, emiprobDist, emiprobVars, nEmis, tprob, 
        iprob);
    em_t *em = setupEM(hmm, emiprobDist, updatetrans, updateemis);
    fwbk_t *fwbk;

    /* number of sequences (chromosomes) and total sequence length */
    int n_seq = Rf_nrows(emi)/ hmm->n_emis;
    int total_seq_length = 0;
    for(int i=0;i<n_seq;i++)
        total_seq_length += Rf_nrows(VECTOR_ELT(emi, i));

/*******************************************************************************
 *
 *  Baum-Welch EM algorithm.
 *
 *  (1) Initialize counts of sufficient stats.
 *  (2) Foreach training sequence (chromosomes) ...
 *    (2a) Run forward/backward algorithms.
 *    (2b) Update sufficient statistics and transition/emission variables.
 *    (2c) Clean up ALL memory that is not absolutely necessary.  Otherwise 
 *    this thing could blow up! :)
 *  (3) Update variables transition and emission variables with new sufficient 
 *      stats.
 *  (4) Calculate the new log-likelihood of the model.
 *  (5) Final cleanup!
 *
 *****************************************************************************/
   if(verb) Rprintf("\nStarting Baum-Welch Algorithm.\n");
      do {
         Q = 0;     // (4)

         /* (1) Initialize counts of sufficient stats. */
         for(int state=0;state<hmm->n_states;state++)  {
            if(em->updateTrans[state]) em->TransSS[state] = 
                    (em->AllocTssFunc[state])(hmm->n_states, n_seq);

           for(int emis=0;emis<hmm->n_emis;emis++)
            if(em->updateEmis[state+hmm->n_states*emis])  
                em->ss[state+hmm->n_states*emis] = 
                    (em->sstats_alloc[state+hmm->n_states*emis])
                        (total_seq_length);
         }

         /* (2) Foreach training sequence (chromosomes) ... */
         for(int seq=0;seq<n_seq;seq++) { // (2)
            // Make the *data variable from the R list type...
            int maxT = Rf_nrows(VECTOR_ELT(emi, seq));
            double **data = (double**)R_alloc(hmm->n_emis, sizeof(double*));
            for(int s=0;s<hmm->n_emis;s++) {
                data[s]   = REAL(VECTOR_ELT(emi, seq+s*n_seq)); 
                //Correctly indexed?!
                assert(Rf_nrows(VECTOR_ELT(emi, seq+s*n_seq))==maxT); 
                // Make sure!
            }
            
            /*  (2a) Run forward/backward algorithms. */
            fwbk= fwbk_alloc(data /* Emission data */, 
                maxT /*(size of *data)*/, hmm);
            forward(fwbk);  // (2a)
            backward(fwbk); // (2a)

            if(verb) Rprintf("Forward prob: %f   Backward prob: %f", 
                fwbk[0].log_px, fwbk[0].bk_log_px);

            /* (2b) Update sufficient statistics and 
             * transition/emission variables. */
            for(int state=0;state<hmm->n_states;state++) {// (2b)
              if(em->updateTrans[state]) (em->UpdateTssFunc[state])(state, 
                seq, em->TransSS[state], hmm->log_eProb, fwbk[0]);

              for(int emis=0;emis<hmm->n_emis;emis++)
                if(em->updateEmis[state+hmm->n_states*emis]) 
                    (em->sstats_emis[state+hmm->n_states*emis])
                    (state, emis, em->ss[state+hmm->n_states*emis], fwbk[0]);
            }
            
            Q+= fwbk[0].log_px;

            /* (2c) Cleanup! */
            fwbk_free(fwbk); // (2c)
         }

         /* (3) Update variables transition and emission variables with 
          * new sufficient stats. */
         if(verb) Rprintf("-- Updating emissions paremeters.\n");
         for(int state=0;state<hmm->n_states;state++) { // (3)
          if(em->updateTrans[state]) (em->UpdateTPFunc[state])(state, n_seq, 
            em->TransSS[state], hmm);

          for(int emis=0;emis<hmm->n_emis;emis++) 
            // DOES THIS ACTUALLY UPDATE SECOND SETS OF EMISSIONS?!
            if(em->updateEmis[state+hmm->n_states*emis]) 
                (em->update_emis[state+hmm->n_states*emis])
                    (state, em->ss[state+hmm->n_states*emis], hmm);
         }
        
         /* (4) Compare change in log-likelihood */
         Qdiff = Q - Qprev; // (4)
         if(verb) Rprintf("-- Likelihood in current, previous (difference) \
            step: %f; %f (%f).\n", Q, Qprev, Qdiff);
         Qprev = Q;     // (4)

         /* (5) Cleanup! */
         for(int state=0;state<hmm->n_states;state++) { // (5) 
            if(em->updateTrans[state])  
                (em->FreeTssFunc[state])(em->TransSS[state]);

          for(int emis=0;emis<hmm->n_emis;emis++)
            if(em->updateEmis[state+hmm->n_states*emis])    
                (em->free_emis_s[state+hmm->n_states*emis])
                    (em->ss[state+hmm->n_states*emis]);  
         }

      } while(Qdiff > T);
   if(verb) Rprintf("EM Converged!  Final log likelihood: %f\n\n", Q);

/*************************************************************
 *  Return data in a list variable to the R enviroment.
 *************************************************************/
    if(verb) Rprintf("Returning to R Enivorment :)\n");
    SEXP ListObject = getEMReturnRTypes(hmm, n_seq, emiprobVars, tprob, emi, 
        output);
    
    return(ListObject);
}

