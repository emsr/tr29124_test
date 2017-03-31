/* cnlnewt.cpp -- Complex Newton solver for nonlinear equations.
 *
 * (C)2001, C. Bond. All rights reserved.
 */
#include <math.h>
#include <complex.h>

void cged(complex<double> **a,complex<double> *b,
    complex<double> *x,int n);

/*
 * Find solution vector for 'f(x) = 0' given initial estimate
 * of x[] and pointer to function which returns residuals fv[].
 *
 * This is a 'bare bones' Newton solver for a set of nonlinear
 * equations. A call to this routine executes a single interation
 * of the Newton method.
 *
 * The Jacobian matrix 'jac' is evaluated numerically using a
 * finite difference approximation based on the user provided
 * 'eps'. Storage for working arrays 'jac' and 'p' is provided
 * by the calling routine, so that allocation overhead needs to
 * be incurred only once in an iterative application.
 */
void cnlnewt(void (*f)(complex<double> *x,complex<double> *fv,int n),
    complex<double> *x,complex<double> *fv,complex<double> **jac,
    complex<double> *p,int n,double eps)
{
    complex<double> tmp,delta;
    int i,j;

    f(x,fv,n);                  // get residuals for current value of 'x'

    for (i=0;i<n;i++) {
        tmp = x[i];
        if (abs(tmp) > 1.0)
            delta = eps*tmp;
        else
            delta = complex<double>(eps,eps);
        x[i] = tmp+delta;       // bump this element
        f(x,p,n);
        x[i] = tmp;             // restore original value
        for (j=0;j<n;j++) {
            jac[j][i] = (p[j]-fv[j])/delta;
        }
    }
    for (i=0;i<n;i++) {
        tmp = 0.0;
        for (j=0;j<n;j++) {
            tmp += jac[i][j]*x[j];
        }
        p[i] = tmp - fv[i];
    }
    cged(jac,p,x,n);

}
/* This version executes several loops of the iteration specified
 * by the user provided variable 'iter'. All auxiliary memory
 * is handled here so the user does not have to manage memory
 * unnecessarily.
 */
void cnlnewt2(void (*f)(complex<double> *x,complex<double> *fv,int n),
    complex<double> *x,complex<double> *fv,int n,double eps,int iter)
{
    complex<double> tmp,delta,**jac,*p;
    int i,j,k;

    p = new complex<double> [n];
    jac = new complex<double> *[n];
    for (i=0;i<n;i++) {
        jac[i] = new complex<double> [n];
    }

    for (k=0;k<iter;k++) {
        f(x,fv,n);                  // get residuals for current value of 'x'

// Compute Jacobian matrix
        for (i=0;i<n;i++) {
            tmp = x[i];
            if (abs(tmp) > 1.0)
                delta = eps*tmp;
            else
                delta = complex<double>(eps,eps);
            x[i] = tmp+delta;       // bump this element
            f(x,p,n);
            x[i] = tmp;             // restore original value
            for (j=0;j<n;j++) {
                jac[j][i] = (p[j]-fv[j])/delta;
            }
        }
// Update residuals
        for (i=0;i<n;i++) {
            tmp = 0.0;
            for (j=0;j<n;j++) {
                tmp += jac[i][j]*x[j];
            }
            p[i] = tmp - fv[i];
        }
// Update solution vector
        cged(jac,p,x,n);
    }
// Delete memory
    for (i=0;i<n;i++) {
        delete [] jac[i];
    }
    delete [] jac;
    delete [] p;
}

/* This version executes up to 'maxiter' loops, exiting if the
 * solution vector converges.
 */
void cnlnewt3(void (*f)(complex<double> *x,complex<double> *fv,int n),
    complex<double> *x,complex<double> *fv,int n,double eps,int *maxiter)
{
    complex<double> tmp,delta,**jac,*p,*x0;
    double dtmp;
    int i,j,k;

    x0 = new complex<double> [n];
    p = new complex<double> [n];
    jac = new complex<double> *[n];
    for (i=0;i<n;i++) {
        jac[i] = new complex<double> [n];
    }

    for (k=0;k<*maxiter;k++) {
        f(x,fv,n);                  // get residuals for current value of 'x'

// Compute Jacobian matrix
        for (i=0;i<n;i++) {
            tmp = x[i];
            if (abs(tmp) > 1.0)
                delta = eps*tmp;
            else
                delta = complex<double>(eps,eps);
            x[i] = tmp+delta;       // bump this element
            f(x,p,n);
            x[i] = tmp;             // restore original value
            for (j=0;j<n;j++) {
                jac[j][i] = (p[j]-fv[j])/delta;
            }
        }
// Update residuals
        for (i=0;i<n;i++) {
            tmp = 0.0;
            for (j=0;j<n;j++) {
                tmp += jac[i][j]*x[j];
            }
            p[i] = tmp - fv[i];
        }
// Update solution vector
        cged(jac,p,x0,n);
        dtmp = 0.0;
        for (i=0;i<n;i++) {
            dtmp += abs(x0[i]-x[i]);
            x[i] = x0[i];
        }
        if (dtmp < 1e-4) break;
    }
    *maxiter = k;
// Delete memory
    for (i=0;i<n;i++) {
        delete [] jac[i];
    }
    delete [] jac;
    delete [] p;
    delete [] x0;
}
