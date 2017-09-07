#ifndef newtH
#define newtH
#include <complex>

void nlnewt(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,double **jac,double *p,int n,double eps);
void nlnewt2(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,int n,double eps,int iter);
void nlnewt3(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,int n,double eps,int *maxiter);
void cnlnewt(void (*f)(std::complex<double> *x,std::complex<double> *fv,int n),
    std::complex<double> *x,std::complex<double> *fv,std::complex<double> **jac,
    std::complex<double> *p,int n,double eps);
void cnlnewt2(void (*f)(std::complex<double> *x,std::complex<double> *fv,int n),
    std::complex<double> *x,std::complex<double> *fv,int n,double eps,int iter);
void cnlnewt3(void (*f)(std::complex<double> *x,std::complex<double> *fv,int n),
    std::complex<double> *x,std::complex<double> *fv,int n,double eps,int *maxiter);

#endif
