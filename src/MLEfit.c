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
 *  Source code written for GRO-seq package by Charles Danko.
 *
 *  Maximum likelihood fit for several distributions
 *
 *  2009-22-11 Maximum likelihood fit for gamma, and R-wrapper for testing.
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

#include "UsefulValues.h"
#include "hmmHeader.h" // Required for expSum and expDif functions.

/*******************************************************************************
 *
 *  diGamma and triGamma functions are needed to calculate MLE of gamma 
 *      distribution.
 *  Numerical approximations are taken from Choi and Wette, 1969.
 *
 *  Commented: 12-3-09; Conflict with R functions of the same purpose.
 *
 ******************************************************************************/
/*
static inline double digamma(double k) {
    if(k < 8) {
        return(digamma(k+1)-1/k);
    }
    else {
        return(log(k)-(1+(1- (0.1-1/(21*k*k)) /(k*k))/(6*k))/(2*k));
    }
}

static inline double trigamma(double k) {
    if(k < 8) {
        return(trigamma(k+1)+1/(k*k));
    }
    else {
        return((1+(1+(1-(1/5-1/(7*k*k))/(k*k))/(3*k))/(2*k))/k);
    }
}
*/

/*******************************************************************************
 *
 *  MLEGamma -- Returns (by reference) maximum likelihood estimates for 
 *  paremeters of the gamma distribution, shape (or k), and scale (or theta).
 *
 *****************************************************************************/
extern int MLEGamma(double N, double SumXi, double SumLogXi, double *shape, 
    double *scale) {
    double shapeBound= VERY_LARGE_DOUBLE_VALUE;
    int maxIterations= 10000;   

    int retVal=0;
    Rprintf("[MLEGamma] SumXi=%f; SumLogXi=%f; N=%f\n", SumXi, SumLogXi, N);

    // Get initial value of shape (k).
    double s = log(SumXi/N)-SumLogXi/N;
    shape[0] = (3-s+sqrt((s-3)*(s-3)+24*s))/(12*s);

    Rprintf("[MLEGamma] s=%f; shape=%f\n", s, shape[0]);

    // Refine shape (k) using Newton's method.
    // int counter = 0; Unused
    double kPrime;
    for( int counter=0; counter < maxIterations; counter++) {
        kPrime= shape[0]-
            (log(shape[0])-digamma(shape[0])-s)/
                ((1/shape[0])-trigamma(shape[0]));
        if((kPrime- shape[0]) < 0.00000001 && (shape[0] - kPrime) <0.00000001){
            shape[0] = kPrime;
            break;
        }

        shape[0]=kPrime;
    }

    // If overflow, use an approximation, and spit out a warning.
    if(isnan(shape[0])) {
        retVal = -1;
        shape[0] = (3-s+sqrt((s-3)*(s-3)+24*s))/(12*s);
        Rprintf("WARNING!! [MLEGamma] NaN returned from Newton's method.  \
            Approximate value returned.\n");
    }

    // (shape > shapeBound) Here, we are bounding shape for numerical stability.
    // (s < 0 && shape < 0) has been observed to happen due to machine rounding 
    // error, if log(SumXi/N)~=SumLogXi/N.
    if(isinf(shape[0]) || (shape[0] > shapeBound) || (s < 0 && shape[0] < 0)){ 
        retVal = -1;
        shape[0] = shapeBound;
        Rprintf("WARNING!! [MLEGamma] Numerical instabillity detected, \
            or shape outside of bounds.\n");
    }

    // Get scale in terms of k.
    scale[0] = (SumXi/(shape[0]*N));

    // This was observed to happen in a case where the first state captured the
    // entire sequence, and N for the second and thrid state --> 0.
    // By setting shape to 0, we are forcing any window that is not 0 to 
    // -Inf prob. 
    // BE CAREFUL!  THIS works ONLY in cases where all windows are set to >=1.
    if(N==0) {
        retVal = -2;
        shape[0]=0; scale[0]=1;
        Rprintf("SERIOUS WARNING!! [MLEGamma] N found to equal 0.  \
            Shape set to 0, scale set to 1.\n");
    }

    Rprintf("[MLEGamma] shape=%f; scale=%f\n", shape[0], scale[0]);
    return(retVal);
}

SEXP RgammaMLE(SEXP n, SEXP sumxi, SEXP sumlogxi) {
    double N = REAL(n)[0];
    double SumXi = REAL(sumxi)[0];
    double SumLogXi = REAL(sumlogxi)[0];

    // Create return value.
    SEXP returnList, returnNames, shape, scale;
    PROTECT(returnList  = allocVector(VECSXP, 2));
    PROTECT(returnNames = NEW_CHARACTER(2));
    SET_VECTOR_ELT(returnList, 0, shape=allocVector(REALSXP, 1));
    SET_VECTOR_ELT(returnList, 1, scale=allocVector(REALSXP, 1));
    SET_STRING_ELT(returnNames, 0, mkChar("shape"));
    SET_STRING_ELT(returnNames, 1, mkChar("scale"));
        setAttrib(returnList, R_NamesSymbol, returnNames);

    MLEGamma(N, SumXi, SumLogXi, REAL(shape), REAL(scale));

    UNPROTECT(2);
    return(returnList);
}

/*******************************************************************************
 *
 *  MLEGammaSHAPE1 -- Returns (by reference) maximum likelihood estimates for 
 *      the hape paremeter of a gamma distribution, where the scale (theta) is 
 *      constrained to 1.
 *
 *  This approximates a poisson where the mean equals the variance.
 *  Scale is set to 1.
 *
 *****************************************************************************/
extern int MLEGamma_SCALE1(double N, double SumXi, double SumLogXi, 
    double *shape, double *scale) {
    int maxIterations= 10000;

    // Get initial value of shape (k).
    double s = SumLogXi/N;
    shape[0] = SumXi/N;

    // Refine shape (k) using Newton's method.
    // int counter = 0; Unused
    double kPrime;
    for( int counter=0; counter < maxIterations; counter++) {
        kPrime=shape[0]-((digamma(shape[0])-s)/trigamma(shape[0]));
        if((kPrime - shape[0]) < 0.00000001 && (shape[0] - kPrime)<0.00000001){
            shape[0] = kPrime;
            break;
        }

        shape[0]=kPrime;
    }

    scale[0]=1;
    return(0);
}

/*******************************************************************************
 *
 *  MLEGamma_SHAPEeq1overSCALE -- Returns (by reference) maximum likelihood 
 *  estimates for the shape paremeter of a gamma distribution, where the scale 
 *  (theta) is constrained to be 1/shape (k).
 *
 ******************************************************************************/
extern int MLEGamma_SHAPEeq1overSCALE(double N, double SumXi, double SumLogXi, 
    double SumXiSq, double *shape, double *scale) {
        int maxIterations= 10000;

    // Get initial value of shape (k).
    double s = (SumXi/N)-(SumLogXi/N);
    shape[0] = ((SumXi/N)*(SumXi/N))/( (SumXiSq/N)-((SumXi/N)*(SumXi/N)) );
    Rprintf("[MLEGamma_SHAPEeq1overSCALE] SumXi=%f; SumLogXi=%f; \
        SumXiSq=%f; N=%f\n", SumXi, SumLogXi, SumXiSq, N);

    // Refine shape (k) using Newton's method.
    // int counter = 0; // Unused
    double kPrime;
    for( int counter=0; counter < maxIterations; counter++) {
        Rprintf("[MLEGamma_SHAPEeq1overSCALE] shape: %f\n", shape[0]);
        kPrime=shape[0]-
            ((digamma(shape[0])+log(1/shape[0])+
                (shape[0]*shape[0])+s)/(trigamma(shape[0])+3*shape[0]));
        if((kPrime - shape[0]) < 0.00000001 && (shape[0] - kPrime)<0.00000001){
            shape[0] = kPrime;
            break;
        }

        shape[0]=kPrime;
    }

    scale[0]=1/shape[0];
    return(0);
}

/*******************************************************************************
 *
 *  normal_exp_optimfn -- Returns the evaluation of the componenet of the 
 *      log-likelihood function for the Normal+Exponential component of the 
 *      Rate HMM.  Used in conjunction with a numerical minimization method 
 *      (currently trying conjugate gradients).
 *
 *      Proof is at the end of the April 2011 - (?) notebook.
 *
 *      As in "Writing R Extensions" Section 6.8: 
 *          typedef double optimfn(int n, double *par, void *ex)
 *              where,  n--> number of args in par.
 *                      par--> parameters.
 *                     ex--> auxiliary information, in this case doubles 
 *                          representing emissions and weights
 *                      encoded as: (1) n, {(2) xi, and (3) log(wi)}_n.
 *
 ******************************************************************************/
extern double normal_exp_optimfn(int n, double *par, void *ex) {
    double *data = (double*)ex; 
    int n_x = (int)(data[0]);
    int maxn = 2*n_x+1, cnt=0;
    double value=0;
    for(int i=1;i<maxn;i+=2) {
        value += data[i+1]*NORMAL_EXP(data[i], par, 4); 
        // cnt is the same as (i-1)/2. Avoid division by keeping counter.
        cnt++;
    }  //   Previously -- worked in log-space.  
      //        double value = expSum(sumCpt, n_x);
    
    return((-1)*value); // Since cgmin minimizes the function, 
                        // have to invert to ID the maximum ... 
}

/*******************************************************************************
 *
 *  normal_exp_optimgr -- Returns (in *gr) the evaluation of the componenet of 
 *          the log-likelihood function for the Normal+Exponential component of 
 *          the Rate HMM.  Used in conjunction with a numerical minimization 
 *          method (currently trying conjugate gradients).
 *
 *          Components of gr are: alpha, mu, sigma, lambda.
 *
 *****************************************************************************/
extern void normal_exp_optimgr(int n, double *par, double *gr, void *ex) {
    double *data = (double*)ex; 
    int n_x = (int)(data[0]);
    int maxn = 2*n_x+1;//, cnt=0;
    double xi, wi, lwi, ximm, N, E, D;

    for(int j=0;j<n;j++)
        gr[j] = 0;

    // int c_lP=0, c_lM=0; Unused
    for(int i=1;i<maxn;i+=2){
      xi  = data[i];
      wi  = data[i+1];
      lwi = log(wi);
      ximm = xi - par[1]; // xi-mu.
      
      // Compute N and E and D
      N = (NORMAL(xi, &par[1], 2));
      E = (EXPONENTIAL(xi, &par[3], 1));
      D = (NORMAL_EXP(xi, par, 4));//expSum2(log(par[0])+N,log(1-par[0])+E); 
      // Compute (alpha*N+(1-alpha)*E).  
      // Same as NORMAL_EXP(), but don't have to recompute N and E.

      gr[0]+= (N>E?1:-1)*exp(lwi+expDif(N,E)-D); 
      // expDif returns the absolute value... but N>0 and E>0, 
      // so abs shoudl be the magnitude.  N>E calculates the sign.
      gr[1]+= ximm*exp(lwi+N-D);
      gr[2]+= ((ximm*ximm)/(par[2]*par[2])-1)*exp(lwi+N-D);
      
    }
    gr[0]= (-1)*gr[0];//exp(expSum(alpha,n_x));
    gr[1]= (-1)*gr[1]*par[0]/par[2]/par[2];
    //exp(par[0]+expSum(mu,n_x)-par[2]-par[2]);
    gr[2]= (-1)*gr[2]*par[0]/par[2];//exp(par[0]+expSum(sigma,n_x)-par[2]);

    // Lots of trouble fitting lambda.  Just use finite difference for 
    // that part...
    double h=0.01, *parl=(double*)calloc(n, sizeof(double*)), 
        *parr=(double*)calloc(n, sizeof(double*));
     for(int j=0;j<n;j++) {
       parl[j] = parr[j] = par[j];
     }
     parl[3] = par[3]-h;
     parr[3] = par[3]+h;
     gr[3]= (normal_exp_optimfn(4, parr, ex)-
        normal_exp_optimfn(4, parl, ex))/(2*h); 
        // DON'T Multiply by -1 for finitDif method.

    return;
}

extern void normal_exp_optimgr_fn_diff_approx(int n, double *par, 
    double *gr, void *ex) {
    double h=0.01, *parl=(double*)calloc(n, sizeof(double*)), 
        *parr=(double*)calloc(n, sizeof(double*));

    for(int i=0;i<n;i++){
     for(int j=0;j<n;j++) {
       parl[j] = parr[j] = par[j];
     }
     parl[i] = par[i]-h;
     parr[i] = par[i]+h;
     gr[i]= (normal_exp_optimfn(4, parr, ex)-
        normal_exp_optimfn(4, parl, ex))/(2*h); 
        // DON'T Multiply by -1 for finitDif method.
    }
    
    return;
}

/******************************************************************************
 *
 *  normal_exp_optcg -- Returns parameters of the (alpha)*normal+(1-alpha)
 *      *exponential distribution optimized using conjugate gradient method.
 *
 *****************************************************************************/ 
SEXP RNormExpMLE(SEXP xi, SEXP wi, SEXP init_guess, SEXP TOL, SEXP maxit) {
    int N = Rf_nrows(xi);
    double *XI = REAL(xi);
    double *WI = REAL(wi);
    double tol= REAL(TOL)[0];
    int MAXIT= INTEGER(maxit)[0];
    double* guess = REAL(init_guess);
    
    // Prepare ex.
    double *ex= (double*)calloc(2*N+1, sizeof(double));
    ex[0] = (double)N;
    for(int cnt=1;cnt<(2*N+1);cnt+=2) {
      ex[cnt]   = XI[(cnt-1)/2];
      ex[cnt+1] = WI[(cnt-1)/2];
    }

    // Create return value.
    SEXP returnList, returnNames, PAR, FMIN, FAIL, FNCOUNT, GRCOUNT;
    PROTECT(returnList  = allocVector(VECSXP, 5));
    PROTECT(returnNames = NEW_CHARACTER(5));
    SET_VECTOR_ELT(returnList, 0, PAR=allocVector(REALSXP, 4));
    SET_VECTOR_ELT(returnList, 1, FMIN=allocVector(REALSXP, 1));
    SET_VECTOR_ELT(returnList, 2, FNCOUNT=allocVector(INTSXP, 1));
    SET_VECTOR_ELT(returnList, 3, GRCOUNT=allocVector(INTSXP, 1));
    SET_VECTOR_ELT(returnList, 4, FAIL=allocVector(INTSXP, 1));
    SET_STRING_ELT(returnNames, 0, mkChar("parameters"));
    SET_STRING_ELT(returnNames, 1, mkChar("minimum.energy"));
    SET_STRING_ELT(returnNames, 2, mkChar("num_function_calls"));
    SET_STRING_ELT(returnNames, 3, mkChar("num_gradient_calls"));
    SET_STRING_ELT(returnNames, 4, mkChar("fail_param"));
    setAttrib(returnList, R_NamesSymbol, returnNames);

    // Convert return values to C compatable format.
    double *par= REAL(PAR);
    double *fmin= REAL(FMIN);
    int *fncount= INTEGER(FNCOUNT);
    int *grcount= INTEGER(GRCOUNT);
    int *fail= INTEGER(FAIL);

    // Do the minimization...
    cgmin(4, guess, par, fmin, normal_exp_optimfn, normal_exp_optimgr, fail, 
        tol, tol, ex, 1, 0, fncount, grcount, MAXIT);

    UNPROTECT(2);
    return(returnList);
}


