/*
 *  glm.h
 *  
 *
 *  Created by Wei Sun on 5/11/10.
 *  Copyright 2010 UNC. All rights reserved.
 *
 */
#ifndef _GLM_H_
#define _GLM_H_

#define BINOMIAL  1
#define POISSON   2
#define GAUSSIAN  3
#define GAMMA     4
#define NB        5
#define ZINB      6

/* Link */

#define LOGIT     1
#define LOG       2
#define IDENTITY  3
#define INVERSE   4

/* GLM definition functions */

double  varfun(int, double, double);
int     muvalid(int, double);
void    initialize(int, double*, double*, int, double*);

/*  log likelihood of Poisson */
double loglik_Poisson(int N, double* mu, double* y);

/* log likelihood of negative binomial */

double loglik_NB(int N, double phi, double* mu, double* y);


/* score and infor for solving MLE of phi */

void score_info(int N, double theta, double* mu, double *y, 
                double* score, double* info);

/* MLE or MM of phi */

int phi_ml(double* y, double* mu, int N,  
           int limit, double eps, double* phi, int initPhi, int trace);

int phi_mm(double* y, double* mu, double* epsR, double* phi, int* dims);

#endif
