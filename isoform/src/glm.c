/*
 *  glm.c
 *  
 *
 *  Created by Wei Sun on 5/11/10.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include "glm.h"

/**********************************************************************
 *
 * utilit functions for glm_fit
 *
 * Originally from R/family and glm_test.cpp of R/CNVTools 
 *
 * Codes invovled strata are removed
 * 
 **********************************************************************/

/* 
 
 Variance function
 
 family:
 1    Binomial
 2    Poisson
 3    Gaussian
 4    Gamma
 5    Negagtive Binomial
 */

double varfun(int family, double mu, double phi){
  switch (family) {
    case 1: return((mu*(1.0-mu)));  /* Binomial */
    case 2: return(mu);             /* Poisson */
    case 3: return(1.0);            /* Gaussian */
    case 4: return(mu*mu);          /* Gamma */
    case 5: return(mu + mu*mu*phi); /* Negative Binomial */
    default: return(0.0);
  }
}

/* Valid values for fitted value, mu. 
 
 If, during iteration, an invalid value is returned, the case is omitted 
 
 */

int muvalid(int family, double mu) {
  double minb = 0.0001, maxb = 0.9999, minp = 0.0001;
  double gammaMin = 0.001;
  switch (family) {
    case 1: return(mu>minb && mu<maxb);    /* Binomial */
    case 2: return(mu>minp);               /* Poisson  */
    case 4: return(mu>gammaMin);           /* Gamma    */
    case 5: return(mu>minp);               /* Negative Binomial */
    default: return(1);                    /* Gaussian */
  }
}

/**********************************************************************
 *
 * utility functions for glm_fit
 *
 * following the implementation as R/stats/family.R 
 * 
 **********************************************************************/

void initialize(int family, double* y, double* mu, int N, double* nTotal_binom) 
{
  int i;
  
  if (family==BINOMIAL) {         /* Binomial */
    for (i=0; i<N; i++) {
      if (y[i] <0 || nTotal_binom[i] < 0) {
        error("negative values not allowed for the Binomial family");
      }
      if (y[i] > nTotal_binom[i]) {
        error("# of success is larger than # of total trials in Binomial family");
      }
      mu[i] = (y[i] + 0.5)/(nTotal_binom[i]+1.0);
    }
  }else if (family==POISSON) {    /* Poisson */
    for (i=0; i<N; i++) {
      if (y[i] < 0) {
        error("negative values not allowed for the Poisson family");
      }      
      mu[i] = y[i] + 0.1;
    }
  }else if (family==GAUSSIAN) {   /* Gaussian*/
    for (i=0; i<N; i++) {
      mu[i] = y[i];
    }
  }else if (family==GAMMA){       /* Gamma */
    for (i=0; i<N; i++) {
      if (y[i] <= 0) {
        error("non-poistive values not allowed for the Gamma family");
      }      
      mu[i] = y[i] + 0.1;
    }
  }else if (family==NB){          /* Negagtive Binomial */
    for (i=0; i<N; i++) {
      if (y[i] < 0) {
        error("negative values not allowed for the Negative Binomial family");
      }else if (y[i] < 0.01) {
        mu[i] = y[i] + 0.1667;
      }else{
        mu[i] = y[i];
      }
    }
  }else {
    error("invaid family");
  }
}

/**********************************************************************
 *
 * log likelihood of Poisson
 *
 **********************************************************************/

double loglik_Poisson(int N, double* mu, double* y){
  int i;
  double yi, mui, logL = 0.0;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    logL += (yi*log(mui) - mui - lgammafn(yi + 1.0));
  }
  
  return(logL);
}

/**********************************************************************
 *
 * log likelihood of negative binomial
 *
 **********************************************************************/

double loglik_NB(int N, double phi, double* mu, double* y){
  int i;
  double logL1, logL, yi, mui;
  double th = 1.0/phi;
  double logL0 = th*log(th) - lgammafn(th);
  
  if (phi < 1e-15) {
    error("loglik_NB: phi is too small\n");
  }
  
  logL = 0.0;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];

    if (yi==0) {
      logL1  = th*log(th) - th*log(th + mui);
    }else {
      logL1  = lgammafn(th + yi) - lgammafn(yi + 1.0) + yi*log(mui) - (th + yi)*log(th + mui);
      logL1 += logL0;
    }

    logL += logL1;
  }
  
  return(logL);
}

/**********************************************************************
 *
 * score_info
 *
 * score and Fisher information, i.e., the first and negative second 
 * derivative of likelihood of phi 
 *
 **********************************************************************/

void score_info(int N, double theta, double* mu, double* y, 
                double* score, double* info)
{
  int i;
  double score1=0.0, info1=0.0;
  double mui, yi, scorei, infoi, thMui;
  
  for (i=0; i<N; i++) {
    yi  = y[i];
    mui = mu[i];
    
    thMui   = theta + mui;
    scorei  = digamma(yi + theta) - digamma(theta) - (theta + yi)/thMui;
    score1 += (scorei - log(thMui) + 1 + log(theta));
    
    infoi   = trigamma(theta) - trigamma(yi + theta) + (mui - yi)/(thMui*thMui);
    info1  += (infoi + 1/thMui - 1/theta);
  }
  
  *score = score1;
  *info  = info1;
}

/**********************************************************************
 *
 * phi_ml
 *
 * MLE of phi (over-dispersion parameter), given mu 
 *
 * Actually we find MLE of 1/phi here and then take inverse
 *
 **********************************************************************/

int phi_ml(double* y, double* mu, int N, int limit, double eps, 
           double* phi, int initPhi, int trace)
{
  double theta0, del, tmp;
  double score=0.0;
  double info=0.0;
  int i, it=0;
  double minTheta = 1e-5;
  double maxTheta = 1.0/minTheta;
  int tryPoisson = 0;
  int tryZINB = 0;
  int fail = 0;
  
  if(initPhi){
    theta0 = 1.0/(*phi);
  }else{
    theta0 = 0.0;
    for (i=0; i<N; i++) {
      tmp = y[i]/mu[i] - 1.0;
      theta0 += tmp*tmp;
    }
    theta0 = (double)N/theta0;
  }
  
  it  = 0;
  del = 1.0;
  
  if(trace > 7) Rprintf("  phi.ml: initial phi = %.2e\n", 1/theta0);
  
  while(it < limit && fabs(del) > eps) {
    score_info(N, theta0, mu, y, &score, &info);
    del     = score/info;
    theta0 += del;
    it     += 1;
    
    if(trace > 7) Rprintf("  phi.ml: iter %d, phi=%.2e, score=%.2e, info=%.2e\n", 
                          it,  1/theta0, score, info);
    
    if (theta0 > maxTheta) {
      theta0 = maxTheta;
      if(trace > 3)
        Rprintf("    phi is truncated at %.2e, no overDispersion?\n", 1/maxTheta);
      
      tryPoisson = 1;
      break;
    }else if (theta0 < 0) {
      
      if(trace > 3)
        Rprintf("    phi will be truncated at 0.0, negative dispersion?\n");
      
      fail = 1;
      break;
      
    }else if (theta0 < minTheta) {
      theta0 = minTheta;
      if(trace > 3)
        Rprintf("    phi is truncated at %.2e, too much overDispersion?\n", 1/minTheta);
      
      tryZINB = 1;
      break;
    }
    
  }
  
  if(it == limit && fabs(del) > 0.01) {
    fail = 1;
    if(trace > 3)
      Rprintf("  phi.ml: iteration limit reached in phi_ml\n");
  }
  
  if(theta0 <= 0){
    *phi = 0.0;
  }else{
    *phi = 1/theta0;
  }
  
  return(tryPoisson + 2*tryZINB + 4*fail);
}

/**********************************************************************
 *
 * phi_mm
 *
 * monment estimator of phi (over-dispersion parameter), given mu 
 *
 * Actually we find MLE of 1/phi here and then take inverse
 *
 **********************************************************************/

int phi_mm(double* y, double* mu, double* epsR, double* phi, int* dims)
{
  int N, dfr, limit, initPhi, trace;
  
  double theta0, del, tmp, mui, dno, eps;
  double minTheta = 1e-5;
  double maxTheta = 1.0/minTheta;
  int tryPoisson = 0;
  int tryZINB = 0;
  int fail = 0;
  int i, it=0;

  N       = dims[0];
  dfr     = dims[1];
  limit   = dims[2];
  initPhi = dims[3];
  trace   = dims[4];
  
  eps     = *epsR;
  
  if(initPhi){
    theta0 = 1.0/(*phi);
  }else{
    theta0 = 0.0;
    for (i=0; i<N; i++) {
      tmp = y[i]/mu[i] - 1.0;
      theta0 += tmp*tmp;
    }
    theta0 = (double)N/theta0;
  }
  
  it  = 0;
  del = 1.0;

  if(trace > 5) Rprintf("  phi.mm: initial phi = %.2e\n", 1/theta0);
    
  while(it < limit && fabs(del) > eps) {
    theta0  = fabs(theta0);
    
    del = dno = 0.0;
    for (i=0; i<N; i++) {
      mui  = mu[i];
      tmp  = (y[i] - mui)*(y[i] - mui);
      del += tmp/(mui + mui*mui/theta0);
      dno += tmp/(mui + theta0)/(mui + theta0);
    }
    
    del    = (del - dfr)/dno;
    theta0 = theta0 - del;
    
    it += 1;
    
    if(trace > 5) Rprintf("  phi.mm: iter %d, phi=%.2e, del=%.2e\n", 
                          it,  1/theta0, del);
    
  }
 
  if (theta0 > maxTheta) {
    theta0 = maxTheta;
    if(trace > 3)
      Rprintf("    phi is truncated at %.2e, no overDispersion?\n", 1/maxTheta);
    
    tryPoisson = 1;
  }else if (theta0 < 0) {
    if(trace > 3)
      Rprintf("    phi will be truncated at 0.0, negative dispersion?\n");
    
    tryPoisson = 1;
  }else if (theta0 < minTheta) {
    theta0 = minTheta;
    if(trace > 3)
      Rprintf("    phi is truncated at %.2e, too much overDispersion?\n", 1/minTheta);
    
    tryZINB = 1;
  }
  
  if(it == limit && fabs(del) > 0.01) {
    fail = 1;
    if(trace > 3)
      Rprintf("  phi.mm: iteration limit reached in phi_ml\n");
  }
  
  if(theta0 <= 0){
    *phi = 0.0;
  }else{
    *phi = 1/theta0;
  }
  
  return(tryPoisson + 2*tryZINB + 4*fail);
}
