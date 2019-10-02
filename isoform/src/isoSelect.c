
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "glm.h"
#include "IAL_nn.h"

/**********************************************************************
 *
 * isoSelect_one
 * isoSelect_one only take one combination of tunning parameters
 * isoSelect calls isoSelect_one search across a number of combinations
 * of tunning parameters 
 *
 
 
 Input:
 
 family       GLM family, eitehr Poisson (2) or NB (5)
 link         Link function (see below)
 N            sample size
 M            # of covariates, i.e., # of columns of X
 y            response, a vector of length N
 offset       offset of invers link function of mean
 X            covariate, a matrix of size N x M
 maxit        maximum number of iterations of IRLS algorithm
 conv         proportional change in weighted sum of squares residuals to
              declare convergence
 init         If true (non-zero), the iteration starts from initial estimates 
              of fitted values (see below). This option has no effect if
              no iteration is required
 
 Output:
 
 fitted       fitted values 
 weights      weights (N-vector)
 beta         regression coeficients
 
 Return
 
 0            convergence
 1            no convergence after maxit iterations
 
 **********************************************************************/

void isoSelect_one(int family, int *dims, double *y, double *muMin,
                   double *offset, double **X, int *succeed,
                   double **wX, double conv, double *fitted,
                   double *weights, double *priorWeights,
                   double *phi, int *trace, int *n2kpIAL, double *BIC, 
                   double lambda1, double tau1, double *wz, 
                   double *b, double *bOld, double *b_new,  
                   double *resIAL, int *w2kpIAL, double *b2kpIAL, 
                   double *Xj2, double *pLm0, double *pLm1)
{
  double epsilon = 1e-15;      /* Singularity threshold */
  double nTotal_binom=0.0;     /* for binomial link, NOT useful here */
  
  int N, M, maxit, maxitIAL, init, useOffset;
  int i, j, k, Nu, p_max, n_lambda, useBIC, usePriorWeights;
  int convged=0, iter=0, useInitValue=0, halfStep, iterHalfStep;
  int maxItHalfStep=3;
  
  if(family!=2 && family!=5){
    Rprintf("family=%d, ", family);
    error("Invalid family! Faimly must be 2 (Poisson) or 5 (NB).\n");
  }
  
  double mu, wi, wsum, vi;
  double xij, *Xj, *wXj, Xjnorm;
  double penalty, logLik=0.0;

  N         = dims[0];
  M         = dims[1];
  maxit     = dims[2];
  maxitIAL  = dims[3];
  init      = dims[4];
  useOffset = dims[5];
  n_lambda  = dims[6];
  p_max     = dims[7];
  useBIC    = dims[8];
  usePriorWeights = dims[9];
  
  if(*trace > 11){
    Rprintf("y\n");
    Rprint_ve(y, 0, 4);
    Rprintf("X\n");
    Rprint_me(X, 0, 1, 0, 4);
  }
  
  if(*trace>3)
    Rprintf("\n  isoSelect_one: family=%d, phi=%f, init=%d\n", 
            family, *phi, init);
    
  /* ----------------------------------------------------------------*
   * by default, initialize mu (fitted) by y itself, with neccesary 
   * modification, e.g., y + 0.1 to avoid log(0) for Poisson family
   * ----------------------------------------------------------------*/
  
  if (init==0) {
    initialize(family, y, fitted, N, &nTotal_binom);
    for (j=0; j<M; j++) { bOld[j] = 0.0; }
  }
  
  /* ----------------------------------------------------------------*
   * Initialize wi (weights) and (standardized) residual 
   *
   * In IRLS, we do a weighted regression z_i ~ X with weight w_i
   *  
   *     z_i = eta_i + (y_i - mu_i)*(d_eta/d_mu)_i
   *   eta_i = Xb = linkfun(mu_i) 
   *   1/w_i = (d_eta/d_mu)^2 v_i
   *
   * where v_i = var(y_i|mu_i,phi)
   *
   * 
   * since we are using identity link function here,
   *      mu = eta
   *     z_i = y_i
   *     w_i = 1/v_i
   *
   * ----------------------------------------------------------------*/
  
  wsum = 0.0;
  
  for (i=0; i<N; i++) {
    
    mu = fitted[i];
    
    if (mu <= 0) { 
      error("invalid mu value\n"); 
    }
    
    vi = varfun(family, mu, *phi);
    wi = 1/vi;
    // later we will look for the number of samples with 
    // non-zeror weight
    if (wi < epsilon) wi = 0.0;

    weights[i] = wi*priorWeights[i];
    
    wsum += wi;
  }
  
  /* ----------------------------------------------------------------*
   * If summation of all weights is too small, stop 
   * ----------------------------------------------------------------*/
  
  if (wsum < epsilon) {
    error("summation of all weights are too small\n"); 
    /*
    if(*trace>2)
      Rprintf("  isoSelect_one: summation of all weights are too small!\n");
    
    maxit    = -1;
    *succeed = 0;
    */
  }else if(*trace > 5){
    Rprintf("\n  isoSelect_one: finish initialization, N=%d, M=%d, family=%d\n", 
            N, M, family);
  }
  
  /* ----------------------------------------------------------------*
   * IRLS algorithm
   *
   * here we adjust the step of Newton's method to make sure
   * the objective function (penalized likelihood) increases
   * ----------------------------------------------------------------*/

  /* set the inital likelihood to be a small number */
  
  convged  = 0;
  iter     = 0;

  while(iter<maxit && !convged) {
    
    if (*trace > 5) {
      Rprintf("\n  isoSelect_one: iteration %d: \n", iter);
    }
        
    /**
     * incoporate the weights to z 
     */
    
    if (useOffset) {
      for (i=0; i<N; i++){
        wz[i] = sqrt(weights[i])*(y[i] - offset[i]);
      }
    }else {
      for (i=0; i<N; i++) {
        wz[i] = sqrt(weights[i])*y[i];
      }
    }    
    
    /**
     * incoporate the weights to X, and
     * calculate the scale of X
     */
    
    for (j=0; j<M; j++) {
      Xj     = X[j];
      wXj    = wX[j];
      Xjnorm = 0.0;
      
      for (i=0; i<N; i++) {
        xij     = sqrt(weights[i])*Xj[i];
        wXj[i]  = xij;
        Xjnorm += xij*xij;
      }
      
      Xj2[j] = Xjnorm;
    }
    
    /* user initial values of b, unless it is the 1st iteration of IRLS
     * and isoSelect_one does not use inital values
     */

    if(init==0 && iter==0){
      useInitValue = 0;
    }else {
      useInitValue = 1;
    }

    IAL_nn(wz, wX, b, b_new, lambda1, tau1, dims, Xj2, conv, resIAL,  
           n2kpIAL, w2kpIAL, b2kpIAL, *trace, useInitValue);
       
    /* for this isoform identificaiton problem, we do not  
     * allow n2kpIAL=0, so exit here and hope another set
     * of tuning parameter can give me better values
     */
    if (*n2kpIAL == 0) {
      *succeed = 0;
      if (*trace > 1) {
        Rprintf("    isoSelect_one: no non-zero covariates is found\n");
      }
      break;
    }
    
    /**
     * double check now, how to adjust the step size of the
     * Newton updates
     */
    
    halfStep     = 1;
    iterHalfStep = 0;
    
    while (halfStep && iterHalfStep < maxItHalfStep) {
      
      /**
       * by default, do not half the step size, 
       * unless we found someting...
       */
      halfStep = 0;

      /* Nu = "N used", the number of samples with non-zero weight */
      Nu = 0;

      /**
       * update the fitted values and weights  
       */
      for (i=0; i<N; i++) {
        mu = 0.0;
        for (k=0; k < *n2kpIAL; k++) {
          j   = w2kpIAL[k];
          mu += X[j][i]*b2kpIAL[k];
        }
        
        if (mu < epsilon) { mu = *muMin; }
        
        if (useOffset) { mu += offset[i]; }
        
        fitted[i] = mu;
                
        vi = varfun(family, mu, *phi);
        /**
         * if vi < epsilon, does not need to update 
         * fitted and weights, just need to half the step size
         */          
        if (vi < epsilon) { 
          halfStep = 1;
          break; 
        }
        
        wi = 1/vi;
        if (wi < epsilon){ wi = 0.0; }else{ Nu ++; }
        weights[i] = wi*priorWeights[i];          
      }
      
      /**
       * if Nu==0, no useful sample, half the step size 
       */          
      
      if (Nu == 0) { halfStep = 1; }
      
      /**
       * if no step size change yet, 
       * check whether penalized likelihood increases
       */          
      
      if(halfStep == 0){
        if (family==2) {
          *pLm1 = loglik_Poisson(N, fitted, y);
        }else {
          *pLm1 = loglik_NB(N, *phi, fitted, y);
        }

        penalty = 0.0;
        for (j=0; j < *n2kpIAL; j++) {
          penalty += log(b2kpIAL[j] + tau1);
        }
        penalty += (M - *n2kpIAL)*log(tau1);
        penalty *= lambda1;

        *pLm1 = *pLm1 - penalty;
        
        if (*pLm1 < *pLm0) {
          if (*trace > 5) {
            Rprintf("    iterHalfStep = %d, (pLm0, pLm1)=(%.2f, %.2f), n2kpIAL=%d\n", 
                    iterHalfStep, *pLm0, *pLm1, *n2kpIAL);
          }      
          
          halfStep = 1; 
        } 
      }
      
      /**
       * if half the step size, 
       * up date b, bOld, w2kpIAL, b2kpIAL, and n2kpIAL
       */          
      if(halfStep == 1){
        iterHalfStep += 1;
        
        k = 0;
        for (j=0; j<M; j++) {
          b[j] = 0.5*b[j] + 0.5*bOld[j];
          
          if(b[j] > 1e-10){ 
            w2kpIAL[k] = j;
            b2kpIAL[k] = b[j];
            k++;
          }
        }
        
        *n2kpIAL = k;
      }

    }
    
    if (*trace > 5) {
      Rprintf("    Finish half stepping, iterHalfStep = %d\n", iterHalfStep);
    }      
    
    if (iterHalfStep >= maxItHalfStep) {
      if (*trace > 5) {
        Rprintf("    isoSelect_one: cannot find good step size in IRLS\n");
      }
      
      /**
       * keep the old coefficients estimates
       */
      
      k = 0;
      for (j=0; j<M; j++) {
        b[j] = bOld[j];
        
        if(b[j] > 1e-10){ 
          w2kpIAL[k] = j;
          b2kpIAL[k] = b[j];
          k++;
        }
      }
      
      *n2kpIAL = k;
      *pLm1 = *pLm0;
      
      //break;
      // cannot break here since otherwise likelihood 
      // is not correct, if this is the first iteration
      // likelihood (logLik) has not be calculated.
      // since b is set as the older value, the algorithm
      // will be claimed to be converged. 
    }
    
    /**
     * update the fitted values  
     */
    for (i=0; i<N; i++) {
      mu = 0.0;
      for (k=0; k < *n2kpIAL; k++) {
        j   = w2kpIAL[k];
        mu += X[j][i]*b2kpIAL[k];
      }
      
      if (mu < epsilon) { mu = *muMin; }
      
      if (useOffset) { mu += offset[i]; }
      
      fitted[i] = mu;
    }
    
    /* record likelihood */

    if (family==2) {
      logLik = loglik_Poisson(N, fitted, y);
    }else {
      logLik = loglik_NB(N, *phi, fitted, y);
    }
    
    penalty = 0.0;
    
    for (j=0; j < *n2kpIAL; j++) {
      penalty += log(b2kpIAL[j] + tau1);
    }
    penalty += (M - *n2kpIAL)*log(tau1);
    penalty *= lambda1;
    
    *pLm1 = logLik - penalty;

    /* finishing up */
    if (*trace > 5) {
      Rprintf("    n2kpIAL=%d, pLm0=%.2f, pLm1=%.2f\n", 
              *n2kpIAL, *pLm0, *pLm1);
    }
    
    //convged  = fabs(*pLm1 - *pLm0)/(fabs(*pLm0) + 0.1) < conv;
    convged  = (*pLm1 - *pLm0)/(fabs(*pLm0) + 0.1) < conv;

    *pLm0 = *pLm1;
    
    for (j=0; j<M; j++) { bOld[j] = b[j]; }
    
    iter++;
  }
  
  if(succeed==0){
    *BIC = 1.0/0.0;
  }else {
    if (useBIC==1 ||(useBIC==0 && M <=N)) {
      //BIC
      *BIC = -2.0*logLik + (*n2kpIAL)*log(N);
    }else{
      // extended BIC
      *BIC = -2.0*logLik + (*n2kpIAL)*(log(N) + log(M));
    }
  }
  
  if(*trace > 2.5){
    Rprintf("    end of isoSelect_one: succeed=%d, n2kpIAL=%d, logLik=%.3e, BIC=%.3e\n", 
            *succeed, *n2kpIAL, logLik, *BIC);
  }
}

/**********************************************************************
 *
 * main function
 *
 **********************************************************************/

void isoSelect(int *dims, double *y, double *muMin, double *offset,
               double *RX, double *convR,  double *scoreTestP, int *trace, 
               double *lambda, double *tau, int *nIter, double *phi2kp, 
               double *priorWeights, int *n2use, int *w2use, int* n2kp, 
               double *b2use, double *Rscore, double *score2use, 
               double *lambda2use, double *tau2use, double *likelihood, 
               int *family)
{
  double MU_Y = 5.0;
  double epsilon = 1e-15;       /* Singularity threshold */

  int N, M, maxit, maxitIAL, init, useOffset, succeed=1;
  int i, j, k1, iter=0, n_lambda, p_max, useBIC, usePriorWeights;
  int fam0=0, n2kpIAL, *w2kpIAL;

  /* dimension parameter for function phi_mm */
  int dims_phi_mm[5];
  
  double conv = *convR, BIC=0.0;
  double del=1.0, pLm0, pLm1, phi0=0.0, lambda1, tau1;
  double Dscore, scoreNum, scoreDen, scorePval, yi, mu;
  double phi[1]; phi[0] = 0.0;
  
  double *wz, **X, **wX;
  double *b, *b_new, *bOld, *Xj2, *resIAL, *b2kpIAL;
  double *fitted, *weights, penalty, *lkhood; 
  
  double sum_mu=0.0, sum_y=0.0;
  
  /* convergence indicator for phi_ml 
   * if cvPhi = 0, NB model is OK. 
   * if cvPhi = 1, suggest we need to use Poisson
   * if cvPhi = 2, suggest we need to use ZINB
   */
  int cvPhi;

  *n2use = 0;
  *score2use = DBL_MAX;
  
  N         = dims[0];
  M         = dims[1];
  maxit     = dims[2];
  maxitIAL  = dims[3];
  init      = dims[4];
  useOffset = dims[5];
  n_lambda  = dims[6];
  p_max     = dims[7];
  useBIC    = dims[8];
  usePriorWeights = dims[9];
  
  /* Each row of X and wX corresponds to a covariate 
   * each column corresponds to an observation
   */
  
  X    = (double **)Calloc(M, double*);
  X[0] = RX;
  for (j=1; j<M; j++) { X[j] = RX + j*N; }
  
  /* wX means weighted X */
  wX    = (double **)Calloc(M, double*);
  wX[0] = (double *) Calloc(M*N, double);
  for(j=1; j<M; j++){ wX[j] = wX[0] + j*N; }

  fitted  = (double *)Calloc(N, double);
  weights = (double *)Calloc(N, double);

  wz = (double *)Calloc(N, double);  // weighted working y
  
  /* variables to be send to function IAL 
   * b is coefficients
   * v is kappa, covaraite-specific weight
   */
  b     = (double *)Calloc(M, double);
  bOld  = (double *)Calloc(M, double);
  b_new = (double *)Calloc(M, double);

  resIAL  = (double *)Calloc(N, double);
  w2kpIAL = (int *)Calloc(M, int);
  b2kpIAL = (double *)Calloc(M, double);
  lkhood  = (double *)Calloc(maxit, double);  // working likelihood

  /* sum square for each marker, i.e., each column of X */
  Xj2 = (double *)Calloc(M, double);
  
  if(*trace > 2.5) {
    Rprintf("\n  isoSelect: N=%d, M=%d, maxit=%d, init=%d, useOffset=%d\n", 
            N, M, maxit, init, useOffset);
  }
  
  for(i=0; i< N; i++){ sum_y += y[i]; }
  
  /* ----------------------------------------------------------------
   * loop across different lambdas and taus
   * ----------------------------------------------------------------*/
  
  for(k1=0; k1 < n_lambda; k1++){
    
    lambda1 = lambda[k1];
    tau1    = tau[k1];
    
    if(*trace > 1){
      Rprintf("\nk1=%d, lambda1=%e, tau1=%e\n", k1, lambda1, tau1);
    }
    
    succeed = 1;
    iter    = 0;

    /* ----------------------------------------------------------
     * Initial fit by Poisson distribution 
     * ----------------------------------------------------------
     */
    
    fam0    = POISSON;
    dims[4] = 0; /* do not use initial values */
    pLm0    = -1.0/epsilon;
    pLm1    = pLm0 + 1.0;
    phi[0]  = 0.0; 
    
    isoSelect_one(fam0, dims, y, muMin, offset, X, &succeed, wX, conv,
                  fitted, weights, priorWeights, phi, trace, &n2kpIAL,
                  &BIC, lambda1, tau1, wz, b, bOld, b_new, resIAL, 
                  w2kpIAL, b2kpIAL, Xj2, &pLm0, &pLm1);
    
    if(!succeed){ 
      if (*trace > 1) {
        Rprintf("\n  fail isoSelect_one\n"); 
      }
      continue; 
    }
    
    /* update fitted value */
    /* test for overdispersion by Dean's Score test */

    scoreNum = 0.0;
    scoreDen = 0.0;
   
    sum_mu = 0.0;
    for (i=0; i<N; i++) {
      /* here mu = fitted[i] since we use identity link function */
      mu = fitted[i];
      yi = y[i];
      sum_mu   += mu;
      scoreNum += (yi - mu)*(yi - mu) - yi;
      scoreDen += mu*mu;
    }
    
    Dscore = scoreNum/sqrt(2.0*scoreDen);
    
    /**
     * double pnorm(double x, double mu, double sigma, 
     *              int lower_tail, int give_log);
     */
    scorePval = pnorm(Dscore, 0.0, 1.0, 0, 0);
    
    if(*trace > 2.5) 
      Rprintf("\n  overdispersion Dean's score = %.2e, p-value = %.2e\n\n", 
              Dscore, scorePval);
    
    /* use Poisson model */
    if(scorePval > *scoreTestP){
      
      if(n2kpIAL > p_max || n2kpIAL == 0){ continue; }
      
      if(sum_mu >= sum_y*MU_Y || sum_mu <= sum_y/MU_Y){
        continue;
      }
      
      Rscore[k1]  = BIC;
      n2kp[k1]    = n2kpIAL;
      
      if (BIC < *score2use) {

        *score2use  = BIC;
        *lambda2use = lambda1;
        *tau2use    = tau1;
        *n2use      = n2kpIAL;
        *nIter      = iter;
        *family     = POISSON;
        *phi2kp     = 0.0;
        
        for (i=0; i < n2kpIAL; i++) {
          w2use[i] = w2kpIAL[i];
          b2use[i] = b2kpIAL[i];
        }
        
        for (i=n2kpIAL; i < p_max; i++) {
          w2use[i] = -9;
          b2use[i] = 0.0;
        }
      }
      
      continue;
    }
    
    /* ----------------------------------------------------------
     * try to use Negative binomial model
     * first estimate dispersion parameter
     * ----------------------------------------------------------
     */
    
    /**
     * calculate phi by MLE, without initial values of phi
     */
          
    cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 0, *trace);
    
    /**
     * if fail phi_ml, try method of moment
     */
    
    if (cvPhi > 3){
      if(*trace > 2.5) 
        Rprintf("\n  fail phi_ml, use phi_mm instead\n");

      dims_phi_mm[0] = N;
      dims_phi_mm[1] = N - n2kpIAL;
      dims_phi_mm[2] = maxit;
      dims_phi_mm[3] = 0;
      dims_phi_mm[4] = *trace;
              
      cvPhi = phi_mm(y, fitted, &conv, phi, dims_phi_mm);
    }
    
    if(cvPhi==0){
      
      if(*trace > 2.5) 
        Rprintf("\n  initial estimate of phi: %e\n", *phi);
      
    }else if (cvPhi==1){
      
      if(*trace > 2.5) 
        Rprintf("\n  choose Poisson family due to small phi:%e\n", *phi);
      
      if(n2kpIAL > p_max || n2kpIAL == 0){ continue; }
      
      if(sum_mu >= sum_y*MU_Y || sum_mu <= sum_y/MU_Y){
        continue;
      }
      
      Rscore[k1]  = BIC;
      n2kp[k1]    = n2kpIAL;
      
      if (BIC < *score2use) {

        *score2use  = BIC;
        *lambda2use = lambda1;
        *tau2use    = tau1;
        *n2use      = n2kpIAL;
        *nIter      = iter;
        *family     = POISSON;
        *phi2kp     = 0.0;

        for (i=0; i < n2kpIAL; i++) {
          w2use[i] = w2kpIAL[i];
          b2use[i] = b2kpIAL[i];
        }
        
        for (i=n2kpIAL; i < p_max; i++) {
          w2use[i] = -9;
          b2use[i] = 0.0;
        }
      }
      
      continue;
      
    }else if(cvPhi==2){
      
      if(*trace > 2.5) 
        Rprintf("\n  the dispersion parameter is too large: phi=%e\n", *phi);
      
      continue;
      
    }else { /* estimation of phi fail to converge */
      if(*trace > 2.5)
        Rprintf("\n  isoSelect: fail to converge in phi_ml: phi=%e\n", *phi);
      
      continue;
    }

    /* ----------------------------------------------------------
     * now, we are sure that we can use Negative binomial model  
     * ----------------------------------------------------------
     */
    
    fam0  = NB;
    
    /**
     * iterative updates
     */
    
    /* del record the differences of phi across iterations */
    
    del  = 1.0; 
    iter = 0;
    
    /* set the inital likelihood to be a small number */
    pLm0 = -1.0/epsilon;
    pLm1 = pLm0 + 1.0;
    
    while (iter < maxit && fabs(pLm0 - pLm1) + fabs(del) > conv) {
      
      if (*trace > 2.5) {
        Rprintf("\n  iteration %d in isoSelect\n", iter);
      }
      
      /* use initial values (1) or not (0) */
      if (iter==0) {
        dims[4] = 0;
      }else {
        dims[4] = 1;
      }
    
      isoSelect_one(fam0, dims, y, muMin, offset, X, &succeed, wX, conv,
                    fitted, weights, priorWeights, phi, trace, &n2kpIAL,
                    &BIC, lambda1, tau1, wz, b, bOld, b_new, resIAL, 
                    w2kpIAL, b2kpIAL, Xj2, &pLm0, &pLm1);
      
      if(!succeed){ 
        if (*trace > 1) {
          Rprintf("\n  fail isoSelect_one\n"); 
        }
        break; 
      }
            
      phi0  = *phi;
      
      /* re-estimate phi */

      cvPhi = phi_ml(y, fitted, N, maxit, conv, phi, 1, *trace);
      
      /**
       * if fail phi_ml, try method of moment
       */
      
      if (cvPhi > 3){
        if(*trace > 2.5) 
          Rprintf("\n  fail phi_ml, use phi_mm instead\n");
        
        dims_phi_mm[0] = N;
        dims_phi_mm[1] = N - n2kpIAL;
        dims_phi_mm[2] = maxit;
        dims_phi_mm[3] = 0;
        dims_phi_mm[4] = *trace;
                  
        cvPhi = phi_mm(y, fitted, &conv, phi, dims_phi_mm);
      }
      
      if(cvPhi==0){
        if(*trace > 2.5) 
          Rprintf("\n  isoSelect: finish phi_ml, phi=%e\n", *phi);
        
      }else if(cvPhi==1){
        if(*trace > 1) 
          Rprintf("  isoSelect: dispersion is too small: phi=%e\n", *phi);
        
        succeed = 0;
        break;
      }else if(cvPhi==2){
        if(*trace > 1) 
          Rprintf("  isoSelect: dispersion is too large: phi=%e\n", *phi);
        
        succeed = 0;
        break;
      }else {
        if(*trace > 1) 
          Rprintf("  isoSelect: fail to converge in phi_ml\n");
        
        succeed = 0;
        break;
      }
      
      del  = phi0 - *phi;
      pLm1 = loglik_NB(N, *phi, fitted, y);
      
      // here the regression coefficients are non-negative
      penalty = 0.0;
      for (j=0; j < n2kpIAL; j++) {
        penalty += log(b2kpIAL[j] + tau1);
      }
      penalty += (M - n2kpIAL)*log(tau1);
      penalty *= lambda1;
  
      if (*trace > 5) {
        Rprintf("\n  pLm1=pLm1 - penalty = %f - %f = %f\n", 
                pLm1, penalty, pLm1 - penalty);
      }
      
      pLm1 = pLm1 - penalty;
      
      if (*trace > 2.5) {
        Rprintf("\n  Phi(%d) = %.2e, pLm1=%f, pLm0=%f, del=%f\n\n", 
                iter, *phi, pLm1, pLm0, del);
      }
      lkhood[iter] = pLm1;
      
      pLm0 = pLm1;
      
      iter++;
    }
    
    if (!succeed) { continue; }
    
    if(iter == maxit) {
      if (*trace > 1) {
        Rprintf("  isoSelect (iter %d): iteration limit reached.\n", iter);
      }
    }
    
    if(n2kpIAL==0 || n2kpIAL > p_max){ 
      if (*trace > 1) {
        Rprintf("  isoSelect (iter %d): skip it due to n2kpIAL=%d\n", iter, n2kpIAL);
      }
      continue; 
    }
    
    sum_mu = 0.0;
    for (i=0; i<N; i++) {
      sum_mu += fitted[i];
    }
    
    
    if(sum_mu >= MU_Y*sum_y){
      if (*trace > 1) {
        Rprintf("  isoSelect (iter %d): skip it since sum_mu(%e) >= sum_y*MU_Y(%e*%e)\n", iter, sum_mu, sum_y, MU_Y);
      }
      continue; 
    }
    
    if(sum_mu <= sum_y/MU_Y){
      if (*trace > 1) {
        Rprintf("  isoSelect (iter %d): skip it since sum_mu(%e) <= sum_y/MU_Y(%e/%e)\n", iter, sum_mu, sum_y, MU_Y);
      }
      continue; 
    }    
    
    Rscore[k1]  = BIC;
    n2kp[k1]    = n2kpIAL;
    
    if(BIC >= *score2use){
      if (*trace > 1) {
        Rprintf("  isoSelect (iter %d): skip it since BIC(%e) < score2use (%e)\n", iter, BIC, *score2use);
      }
      continue; 
    }
        
    if (*trace > 1) {
      Rprintf("  isoSelect (iter %d): done. \n", iter);
      Rprintf("\n  n2kpIAL=%d, BIC=%e\n", n2kpIAL, BIC);
      if (*trace > 2 && n2kpIAL > 0) {
        Rprintf("  w2kpIAL: ");
        Rprint_vi(w2kpIAL, 0, n2kpIAL-1);
        Rprintf("  b2kpIAL: ");
        Rprint_ve(b2kpIAL, 0, n2kpIAL-1);          
      }
    }

    if (BIC < *score2use && sum_mu < sum_y*MU_Y && sum_mu > sum_y/MU_Y) {
      *score2use  = BIC;
      *lambda2use = lambda1;
      *tau2use    = tau1;
      *n2use      = n2kpIAL;
      *nIter      = iter;
      *family     = NB;
      *phi2kp     = *phi;
              
      for (i=0; i < n2kpIAL; i++) {
        w2use[i] = w2kpIAL[i];
        b2use[i] = b2kpIAL[i];
      }
      
      for (i=n2kpIAL; i < p_max; i++) {
        w2use[i] = -9;
        b2use[i] = 0.0;
      }

      for (i=0; i < iter; i++) {
        likelihood[i] = lkhood[i];
      }

      for (i=iter; i < maxit; i++) {
        likelihood[i] = 0.0;
      }
      
    }
  }// end of loop for lambda and tau
   
  
  Free(wz);
  Free(b);
  Free(bOld);
  Free(b_new);
  Free(resIAL);
  Free(w2kpIAL);
  Free(b2kpIAL);
  Free(Xj2);
  Free(lkhood);
  
  Free(X);

  Free(wX[0]);
  Free(wX);
  
  Free(fitted);
  Free(weights);
}
