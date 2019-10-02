/*
 *  IAL.c
 *
 *  Created by Wei Sun on Sep 8 2010.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <R.h>
#include <Rmath.h>
#include "IAL_nn.h"
#include <R_ext/Applic.h>

/**********************************************************************
 * 
 * utility functions
 *
 **********************************************************************/

void Rprint_vi(int* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%d\t", v[i]);
	}
	Rprintf ("%d\n", v[i]);
}

void Rprint_ve(double* v, int nrl, int nrh)
{
	int i;
	for (i = nrl; i < nrh; i++){
		Rprintf ("%.2e\t", v[i]);
	}
	Rprintf ("%.2e\n", v[i]);
}

void Rprint_me(double** m, long nrl, long nrh, long ncl, long nch)
{
	int i, j;
	for (i = nrl; i <= nrh; i++){
    for(j = ncl; j < nch; j++){
      Rprintf ("%.2e\t", m[i][j]);
    }
    Rprintf ("%.2e\n", m[i][j]);
	}
}

static int checkConvergence(double *beta, double *beta_old, double eps, int p)
{
  int j;
  int converged = 1;
  
  for (j=0; j < p; j++)
  {
    if (beta[j]!=0 && beta_old[j]!=0)
    {
      if (fabs((beta[j]-beta_old[j])/beta_old[j]) > eps)
	    {
	      converged = 0;
	      break;
	    }
    }
    else if (beta[j]==0 && beta_old[j]!=0)
    {
      converged = 0;
      break;
    }
    else if (beta[j]!=0 && beta_old[j]==0)
    {
      converged = 0;
      break;
    }
  }
  
  return(converged);
}

/*********************************************************************
 *
 * IAL_nn (nn means non-negative constraints)
 *
 * The Iterative Adaptive Lasso with additionl restrction that 
 * all coef should be non-negative and there is no intercept
 
   y:       response varaible
   X:       covaraite marix, one row for one covariate
   b:       regression coefficient except intercept
   lambda:  tuning parameters lambda
   tau:     tuning parameters tau
   dims:    dimensions
   conv:    a small positive number of to call convergence
   resid:   working residuals, has been initialized by y
   Xj2:     sum square of Xj
   n2kp:    number of covaraites seleted by IAL
 
 *********************************************************************/

void IAL_nn(double* y, double** X, double* b, double* b_new, 
            double lambda1, double tau1, int* dims, 
            double *Xj2, double conv, double* resid, 
            int* n2kp, int* w2kp, double* b2kp, int trace, 
            int useInitValue)
{
  int i, j, k, w, n, p, maxitIAL, converged=0;
  double bj_bar, threshold, bj0, bj, dbj, *Xj;
  
  n = dims[0];
  p = dims[1];
  maxitIAL = dims[3];
    
  if(trace > 7){
    Rprintf("  IAL: n=%d, p=%d, maxitIAL=%d, useInitValue=%d\n", 
            n, p, maxitIAL, useInitValue);
  }

  /**
   * step 1. Initialization of bj
   */
  
  if(! useInitValue){
    for(j=0; j<p; j++){ 
      b[j] = 0.0;
    }
  }
  
  /**
   * step 2 Initialization of resid
   */
    
  for(i=0; i<n; i++){ 
    resid[i] = y[i];
  }
  
  if(useInitValue){
    
    for(j=0; j<p; j++){
      Xj = X[j];
      bj = b[j];
      
      if (bj > 1e-16) {
        for(i=0; i<n; i++){ resid[i] -= Xj[i]*bj; }
      }
    }
    
  }
  
  if(trace > 13){
    Rprintf("y\n");
    Rprint_ve(y, 0, p-1);
    Rprintf("resid\n");
    Rprint_ve(resid, 0, p-1);    
    Rprintf("X\n");
    Rprint_me(X, 0, 1, 0, p-1);
    Rprintf("b\n");
    Rprint_ve(b, 0, p-1);
  }
  
  /* ------------------------------------------------
   * Iterative updating
   * ------------------------------------------------
   */
  
  for(w=0; w<maxitIAL; w++){
    
    /**
     * step 1 Update b[j]
     */
    for(j=0; j<p; j++){
      Xj  = X[j];
      bj0 = b[j];
      
      /* recalculate regression coefficient */
      bj_bar = 0.0;

      for(i=0; i<n; i++){
        bj_bar += Xj[i]*resid[i];
      }
      
      bj_bar = bj_bar/Xj2[j] + bj0;
      
      threshold = lambda1/(bj0 + tau1);
      
      if(bj_bar > threshold){
        bj = bj_bar - threshold;
      }else{
        bj = 0.0;
      }

      dbj = bj - bj0;
      
      /* adjust the residuals based on the new coefficient */
      if (fabs(dbj) > 1e-16) {
        for(i=0; i<n; i++){ resid[i] -= Xj[i]*dbj; }
      }
      
      b_new[j] = bj;
    }
    
    /**
     * check convergence
     */          
	  if (checkConvergence(b_new, b, conv, p))
    {
      converged  = 1;
      break;
    }
    
    if(trace > 11){
      Rprintf("      w=%d, converged=%d\n", w, converged);
    }
    
    for (j=0; j<p; j++) b[j] = b_new[j];

    if(converged){ break; }
    
  }

  /* ------------------------------------------------
   * record results
   * ------------------------------------------------
   */
  
  k = 0;
  
  for(j=0; j < p; j++){
    if(b[j] > 1e-10){ 
      w2kp[k] = j;
      b2kp[k] = b[j];
      k++; 
    }
  }
  
  if (trace > 9) {
    Rprintf("    n2kp=%d\n", k);
    if (k > 0) {
      Rprintf("    w2kp=");
      Rprint_vi(w2kp, 0, k-1);
      Rprintf("    b2kp=");
      Rprint_ve(b2kp, 0, k-1);
    }
  }

  *n2kp = k;
  
}
