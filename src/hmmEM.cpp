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


/********************************************************************************
 *
 *	Baum Welch EM implementation -- Written for the GRO-seq package by 
 *	Charles Danko.
 *
 *	TODO: 
 *	  (1) Look into R support for multi-threading.
 *	  (2) ...
 *
 *	2009-11-23 Imported and revised Forward/backward implementation by Andre Martins.
 *	2009-11-24 Implemented Baum Welch algorthm.
 *
 ********************************************************************************/

using namespace std;

extern "C" {

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

/**********************************************************************************************
 *
 * rBaumWelchEM -- Implementation of the Baum Welch Expectation Maximazation algorithm used to 
 *	estimate the paremeters of an HMM.
 *
 *  nSeq        -- 
 *	nstates	    --	Number of states in the HMM.
 *	emi	        --	List of real vectors representing the emissions over different sequences (e.g. chromosomes).
 *  seqLength   --  1x n ... Vector of sequence lengths.
 *  nEmis       --  Number of emission probabilities for each sequence.
 *	emiprobD    -- 	Vector (1 x nstates) of strings representing the emission probability distirubitons.
 *	emiprobV    --	Matrix (3 x nstates) representing arguments to pdist*.  Unused paremeters can be set to 0.
 *	tprob	    --	Matrix (nstates x nstates) representing transition probabilities between states.
 *	iprob 	    --	Vector (1 x nstates) of initial probabilities.
 *	threshold   --  Real. Threshold for convergence.
 *	updateTrans --	Boolean array.  Rows in the transition parameter matrix that should be updated by EM.
 *	updateEmis  --	Boolean... [TODO: This should be bool* -- which state's params are updated.]
 *  output      --  Integer.  Controls the output list options.
 *				=0|1 -- Standard, data+viterbi.
 *				>1 -- Returns posteriors.
 *				=10-- Returns posteriors & posterior prob. of transition from state 2->3.
 *	verbose     --  Integer. Controls the amount of output text.  
 *				=0 -- quiet, no output; 
 *				=1 -- regular updates;
 *
 *	Questions/Implementation Notes:
 *	-- Can be used to update transition probabilities, emission probabilities, or both using the update*
 *	   variables.
 *	-- Currently no support for using multiple threads on multithread processors.  This may come in
 *	   the future depending on how support is done in R.  It is also possible to implement this support
 *	   on the C-side, but will require additional dependencies.
 *	-- Removes data from each hmm run as quickly as possible to prevent memory from blowing up!
 *	
 *	Updates/Information:
 *	2009-11-24 Wrote this function.
 *
 **********************************************************************************************/
SEXP RBaumWelchEM(/*SEXP nSeq,*/ SEXP nstates, SEXP emi, /*SEXP seqLength*/ SEXP nEmis, SEXP emiprobDist, SEXP emiprobVars, SEXP tprob, SEXP iprob, 
			SEXP threshold, SEXP updatetrans, SEXP updateemis, SEXP output, SEXP verbose) {

/*************************************************************
 *  Set up data structures that are used in EM.
 *************************************************************/

	/* Init logical variables. */
	bool verb = INTEGER(verbose)[0];
	int outopt = INTEGER(output)[0];

	if(verb) Rprintf("Initializing Baum-Welch EM.\n");

	/* booleans indicating whether or not to update by EM */
	int *updateTrans = INTEGER(updatetrans);
	int *updateEmis  = INTEGER(updateemis);

	/* likelihoods */
	double Qprev = -999999999999999.0;
	double Q, Qdiff;
	double T = REAL(threshold)[0];

	/* transfer info from R vars into hmm struct */
	hmm_t *hmm = setupHMM(nstates, emiprobDist, emiprobVars, nEmis, tprob, iprob);
    fwbk_t *fwbk;

	/* number of sequences (chromosomes) and total sequence length */
	int n_seq = /*INTEGER(nSeq)[0];//*/Rf_nrows(emi)/ hmm[0].n_emis;
	int total_seq_length = 0;
	for(int i=0;i<n_seq;i++)
		total_seq_length += /*seq_length[i];//*/Rf_nrows(VECTOR_ELT(emi, i));

	if(verb) 
	  for(int k=0;k<hmm[0].n_states;k++)
		for(int l=0;l<hmm[0].n_states;l++)
			Rprintf("[BWem] k=%d; l=%d; Tprob=%f\n", k, l, hmm[0].log_tProb[k][l]);

	/***************************************************
	* Set up the transition probability distributions.   
	****************************************************/
	alloc_trans_sstats  AllocTssFunc = (alloc_trans_sstats) R_alloc(hmm[0].n_states, sizeof( alloc_trans_sstats ));
	update_trans_SS    UpdateTssFunc = (update_trans_SS)    R_alloc(hmm[0].n_states, sizeof( update_trans_SS ));
	update_trans_Prob   UpdateTPFunc = (update_trans_Prob)  R_alloc(hmm[0].n_states, sizeof( update_trans_Prob ));
	free_trans_sstats    FreeTssFunc = (free_trans_sstats)  R_alloc(hmm[0].n_states, sizeof( free_trans_sstats ));
	void** TransSS = (void**)R_alloc(hmm[0].n_states, sizeof(void*));

	for(int i=0;i<(hmm[0].n_states);i++) {
	  AllocTssFunc[i]   = TransAlloc;
	  UpdateTssFunc[i]  = TransUpdate;
	  UpdateTPFunc[i]   = TransUpdateP;
	  FreeTssFunc[i]    = TransFree;
	}

	/***************************************************
	* Set up the emission probability distributions.   
	****************************************************/
	alloc_emis_sstats sstats_alloc = (alloc_emis_sstats) R_alloc(hmm[0].n_states*hmm[0].n_emis, sizeof( alloc_emis_sstats ));
	update_sstat_func sstats_emis  = (update_sstat_func) R_alloc(hmm[0].n_states*hmm[0].n_emis, sizeof( update_sstat_func ));
	update_emiss_func update_emis  = (update_emiss_func) R_alloc(hmm[0].n_states*hmm[0].n_emis, sizeof( update_emiss_func ));
	free_emis_sstats  free_emis_s  =  (free_emis_sstats) R_alloc(hmm[0].n_states*hmm[0].n_emis, sizeof( free_emis_sstats  ));
	void** ss = (void**)R_alloc(hmm[0].n_states, sizeof(void*));

	for(int i=0;i<(hmm[0].n_states*hmm[0].n_emis);i++) {
		if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "norm") == 0)	{
			sstats_alloc[i] = SSallocNormal;
			sstats_emis[i]  = SStatsNormal;
			update_emis[i]  = UpdateNormal;
			free_emis_s[i]  = SSfreeNormal;
		}
		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "dnorm") == 0)	{
			sstats_alloc[i] = SSallocNormal;
			sstats_emis[i]  = SStatsNormal;
			update_emis[i]  = UpdateNormal;
			free_emis_s[i]  = SSfreeNormal;
		}
		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "gamma") == 0)	{
			sstats_alloc[i] = SSallocGamma;
			sstats_emis[i]  = SStatsGamma;
			update_emis[i]  = UpdateGamma;
			free_emis_s[i]  = SSfreeGamma;
		}
		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "dgamma") == 0)	{
			sstats_alloc[i] = SSallocGamma;
			sstats_emis[i]  = SStatsGamma;
			update_emis[i]  = UpdateGamma;
			free_emis_s[i]  = SSfreeGamma;
		}
		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "normexp") == 0)	{
			sstats_alloc[i] = SSallocNormExp;
			sstats_emis[i]  = SStatsNormExp;
			update_emis[i]  = UpdateNormExp;
			free_emis_s[i]  = SSfreeNormExp;
		}
		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "normexpminus") == 0)	{
			sstats_alloc[i] = SSallocNormExp;
			sstats_emis[i]  = SStatsNormExp;
			update_emis[i]  = UpdateNormExp;
			free_emis_s[i]  = SSfreeNormExp;
		}
/*		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "gamma_scale1") == 0)	{
			sstats_alloc[i] = SSallocGamma;
			sstats_emis[i]  = SStatsGamma;
			update_emis[i]  = UpdateGamma_SCALE1;
			free_emis_s[i]  = SSfreeGamma;
		}
		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "gamma_SHAPEeq1overSCALE") == 0)	{
			sstats_alloc[i] = SSallocGamma;
			sstats_emis[i]  = SStatsGamma;
			update_emis[i]  = UpdateGamma_SHAPEeq1overSCALE;
			free_emis_s[i]  = SSfreeGamma;
		}
		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "gamma_p1") == 0)	{
			sstats_alloc[i] = SSallocGamma;
			sstats_emis[i]  = SStatsGamma_p1;
			update_emis[i]  = UpdateGamma;
			free_emis_s[i]  = SSfreeGamma;
		}*/
		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "pois") == 0)	{
			sstats_alloc[i] = SSallocPoisson;
			sstats_emis[i]  = SStatsPoisson;
			update_emis[i]  = UpdatePoisson;
			free_emis_s[i]  = SSfreePoisson;		
		}
		else if(strcmp(CHAR(STRING_ELT(emiprobDist, i)), "exp") == 0)	{
			sstats_alloc[i] = SSallocExp;
			sstats_emis[i]  = SStatsExp;
			update_emis[i]  = UpdateExp;
			free_emis_s[i]  = SSfreeExp;		
		}
		else error("Distribution ('%s') not recognized!", CHAR(STRING_ELT(emiprobDist, i)));
	}

/**********************************************************************************
 *
 *  Baum-Welch EM algorithm.
 *
 *	(1) Initialize counts of sufficient stats.
 *	(2) Foreach training sequence (chromosomes) ...
 *	  (2a) Run forward/backward algorithms.
 * 	  (2b) Update sufficient statistics and transition/emission variables.
 *	  (2c) Clean up ALL memory that is not absolutely necessary.  Otherwise this thing could blow up! :)
 *	(3) Update variables transition and emission variables with new sufficient stats.
 *	(4) Calculate the new log-likelihood of the model.
 *  (5) Final cleanup!
 *
 **********************************************************************************/
   if(verb) Rprintf("\nStarting Baum-Welch Algorithm.\n");
      do {
         Q = 0;		// (4)

		 /* (1) Initialize counts of sufficient stats. */
		 //Rprintf(" (1) Initialize counts of sufficient stats. \n");
         for(int state=0;state<hmm[0].n_states;state++)  {
			if(updateTrans[state]) TransSS[state] = (AllocTssFunc[state])(hmm[0].n_states, n_seq);

		   for(int emis=0;emis<hmm[0].n_emis;emis++)
            if(updateEmis[state+hmm[0].n_states*emis])  ss[state+hmm[0].n_states*emis] = (sstats_alloc[state+hmm[0].n_states*emis])(total_seq_length);
		 }

		 /* (2) Foreach training sequence (chromosomes) ... */
		 //Rprintf(" (2) Foreach training sequence (chromosomes) ... \n");
         for(int seq=0;seq<n_seq;seq++) { // (2)
            // Make the *data variable from the R list type...
            int maxT = Rf_nrows(VECTOR_ELT(emi, seq));//INTEGER(seqLength)[seq];//Rf_nrows(VECTOR_ELT(emi, seq));
			double **data = (double**)R_alloc(hmm[0].n_emis, sizeof(double*));
			for(int s=0;s<hmm[0].n_emis;s++) {
				data[s]   = REAL(VECTOR_ELT(emi, seq+s*n_seq)); //Correctly indexed?!
				assert(Rf_nrows(VECTOR_ELT(emi, seq+s*n_seq))==maxT); // Make sure!
			}
		
//            double *data = REAL(VECTOR_ELT(emi, seq)); 
			
			/*	(2a) Run forward/backward algorithms. */
			//Rprintf("  (2a) Run forward/backward algorithms... Allocate\n");
            fwbk= fwbk_alloc(data /* Emission data */, maxT /*(size of *data)*/, hmm);
			//Rprintf("  (2b) Run forward.\n");
            forward(fwbk);  // (2a)
			//Rprintf("  (2c) Run backward.\n");
            backward(fwbk); // (2a)

			if(verb) Rprintf("Forward prob: %f   Backward prob: %f", fwbk[0].log_px, fwbk[0].bk_log_px);

			/* (2b) Update sufficient statistics and transition/emission variables. */
            for(int state=0;state<hmm[0].n_states;state++) {// (2b)
              if(updateTrans[state]) (UpdateTssFunc[state])(state, seq, TransSS[state], hmm[0].log_eProb, fwbk[0]);

			  for(int emis=0;emis<hmm[0].n_emis;emis++)
                if(updateEmis[state+hmm[0].n_states*emis]) (sstats_emis[state+hmm[0].n_states*emis])(state, emis, ss[state+hmm[0].n_states*emis], fwbk[0]);
			}
			
            Q+= fwbk[0].log_px;
//            if(verb) Rprintf("-- Likelihood on seq %d= %f;\n", seq, fwbk[0].log_px);

			/* (2c) Cleanup! */
            fwbk_free(fwbk); // (2c)
         }

		 /* (3) Update variables transition and emission variables with new sufficient stats. */
         if(verb) Rprintf("-- Updating emissions paremeters.\n");
         for(int state=0;state<hmm[0].n_states;state++) { // (3)
          if(updateTrans[state]) (UpdateTPFunc[state])(state, n_seq, TransSS[state], hmm);

		  for(int emis=0;emis<hmm[0].n_emis;emis++) // DOES THIS ACTUALLY UPDATE SECOND SETS OF EMISSIONS?!
 		    if(updateEmis[state+hmm[0].n_states*emis]) (update_emis[state+hmm[0].n_states*emis])(state, ss[state+hmm[0].n_states*emis], hmm);
		 }
		
		 /* (4) Compare change in log-likelihood */
         Qdiff = Q - Qprev;	// (4)
         if(verb) Rprintf("-- Likelihood in current, previous (difference) step: %f; %f (%f).\n", Q, Qprev, Qdiff);
         Qprev = Q;		// (4)

		 /* (5) Cleanup! */
         for(int state=0;state<hmm[0].n_states;state++) { // (5) 
            if(updateTrans[state])	(FreeTssFunc[state])(TransSS[state]);

          for(int emis=0;emis<hmm[0].n_emis;emis++)
            if(updateEmis[state+hmm[0].n_states*emis]) 	(free_emis_s[state+hmm[0].n_states*emis])(ss[state+hmm[0].n_states*emis]);	
		 }

      } while(Qdiff > T);
   if(verb) Rprintf("EM Converged!  Final log likelihood: %f\n\n", Q);

/*************************************************************
 *  Return data in a list variable to the R enviroment.
 *************************************************************/
	if(verb) Rprintf("Returning to R Enivorment :)\n");
	// Add final model paremeters to the list.
	int RETURN_SIZE = 3, RETURN_INDX = 0;
	if(outopt>1) RETURN_SIZE++;
	if(outopt==10) RETURN_SIZE++;

	SEXP ListObject, hiddenStatesR, posteriors, postTrans; // ListObject can be named more easily in the R function that calls it.
	PROTECT(ListObject = allocVector(VECSXP, RETURN_SIZE));
	SET_VECTOR_ELT(ListObject, RETURN_INDX++, emiprobVars);
	SET_VECTOR_ELT(ListObject, RETURN_INDX++, tprob);
//	SET_VECTOR_ELT(ListObject, 2, Q);

	// Run Viterbi to get the most likely state paths?!
	SET_VECTOR_ELT(ListObject, RETURN_INDX++, hiddenStatesR=allocVector(VECSXP, n_seq));
	if(outopt > 1) {
//		SEXP posteriors;
		SET_VECTOR_ELT(ListObject, RETURN_INDX++, posteriors=allocVector(VECSXP, n_seq));
	}
	if(outopt == 10) {
		SET_VECTOR_ELT(ListObject, RETURN_INDX++, postTrans=allocVector(VECSXP, n_seq));
	}

	for(int seq=0;seq<n_seq;seq++) {

		// Data in types.
		int maxT = Rf_nrows(VECTOR_ELT(emi, seq));//INTEGER(seqLength)[seq];//Rf_nrows(VECTOR_ELT(emi, seq));
		double **data = (double**)R_alloc(hmm[0].n_emis, sizeof(double*));
		for(int s=0;s<hmm[0].n_emis;s++) {
			data[s]   = REAL(VECTOR_ELT(emi, seq+s*n_seq)); //Correctly indexed?!
			assert(Rf_nrows(VECTOR_ELT(emi, seq+s*n_seq))==maxT); // Make sure!
		}
		
		// Return types.
		SET_VECTOR_ELT(hiddenStatesR, seq, allocVector(INTSXP,maxT));
		int *hiddenstates = INTEGER(VECTOR_ELT(hiddenStatesR, seq));

		viterbi_path(hmm[0], data, maxT, NULL, NULL, hiddenstates);

		// CGD --> 6-6-2011 --> Optionally also returns posteriors (only one sequence currently supported).
		if(outopt > 1) {
			// Set up R data types.  
			// How to set up a matrix?!
			SET_VECTOR_ELT(posteriors, seq, allocMatrix(REALSXP, hmm[0].n_states, maxT)); // REALSXP an actual type?!
			double *Rpost = REAL(VECTOR_ELT(posteriors, seq));
			
			// Run forward/ backward.
			fwbk= fwbk_alloc(data /* Emission data */, maxT /*(size of *data)*/, hmm);
			forward(fwbk);  // (2a)
			backward(fwbk); // (2a)
			Q= fwbk[0].log_px;

			// Foreach state, calculate posteriors.  Construction similar to: http://stat.ethz.ch/R-manual/R-devel/doc/manual/R-exts.html#Attributes
			for(int state=0;state<hmm[0].n_states;state++)
			  for(int position=0;position<maxT;position++)
				Rpost[state + hmm[0].n_states*position] = fwbk[0].forward[position][state]+fwbk[0].backward[position][state]-Q; // Durbin book, eq. 3.14.
				
			if(outopt == 10) { // Return post. prob. of transition between states 2 and 3.
				SET_VECTOR_ELT(postTrans, seq, allocVector(REALSXP, (maxT-1)));
				double *post_trans = REAL(VECTOR_ELT(postTrans, seq));
				for(int position=0;position<(maxT-1);position++) {
				  post_trans[position] = fwbk[0].forward[position][1]+fwbk[0].backward[position+1][2]+hmm[0].log_tProb[1][2] -Q; // +emission // Durbin book, eq. 3.19.
  		 		  for(int emis_count=0;emis_count<hmm[0].n_emis;emis_count++) 
 				    post_trans[position] += (hmm[0].log_eProb[2+hmm[0].n_emis*emis_count])(data[emis_count][position+1], hmm[0].em_args[2+hmm[0].n_emis*emis_count], 4);
				}
			}
			
			// Cleanup...
			fwbk_free(fwbk);
	   }
	}
	
	// Return list object.
	unprotect(1);
	return(ListObject);
}

}
