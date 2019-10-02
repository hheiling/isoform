
void Rprint_vi(int* v, int nrl, int nrh);

void Rprint_ve(double* v, int nrl, int nrh);

void Rprint_me(double** m, long nrl, long nrh, long ncl, long nch);

void IAL_nn(double* y, double** X, double* b, double* b_new, 
            double lambda1, double tau1, int* dims, 
            double *Xj2, double conv, double* resid, 
            int* n2kp, int* w2kp, double* b2kp, int trace, 
            int useInitValue);
