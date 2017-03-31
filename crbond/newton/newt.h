#ifndef newtH
#define newtH
#include <complex.h>

void nlnewt(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,double **jac,double *p,int n,double eps);
void nlnewt2(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,int n,double eps,int iter);
void nlnewt3(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,int n,double eps,int *maxiter);
void cnlnewt(void (*f)(complex<double> *x,complex<double> *fv,int n),
    complex<double> *x,complex<double> *fv,complex<double> **jac,
    complex<double> *p,int n,double eps);
void cnlnewt2(void (*f)(complex<double> *x,complex<double> *fv,int n),
    complex<double> *x,complex<double> *fv,int n,double eps,int iter);
void cnlnewt3(void (*f)(complex<double> *x,complex<double> *fv,int n),
    complex<double> *x,complex<double> *fv,int n,double eps,int *maxiter);

#endif
