/* nlnewt.cpp -- Newton solver for nonlinear equations.
 *
 * (C)2001, C. Bond. All rights reserved.
 */
#include <cmath>
#include <iostream>

int gelimd(double **a,double *b,double *x,int n);

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

void nlnewt(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,double **jac,double *p,int n,double eps)
{
    double tmp,delta;
    int i,j;

    
    f(x,fv,n);                  // get residuals for current value of 'x'

// Compute Jacobian matrix
    for (i=0;i<n;i++) {
        tmp = x[i];
        delta = (tmp > 1.0) ? eps*tmp : eps;
        x[i] = tmp+delta;       // bump this element
        delta = x[i] - tmp;     // try this to reduce error (from Schnabel)
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
    gelimd(jac,p,x,n);
}

/* The following version executes several loops of the iteration
 * specified by the user provided variable 'iter'. All temporary
 * storage (for 'jac' and 'p') are handled here so the user
 * does not have to manage memory unnecessarily.
 */
void nlnewt2(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,int n,double eps,int iter)
{
    double tmp,delta,**jac,*p;
    int i,j,k;

    jac = new double *[n];
    for (i=0;i<n;i++) {
        jac[i] = new double [n];
    }
    p = new double [n];
    
    for (k=0;k<iter;k++) {
        f(x,fv,n);                  // get residuals for current value of 'x'

// Compute Jacobian matrix
        for (i=0;i<n;i++) {
            tmp = x[i];
            delta = (tmp > 1.0) ? eps*tmp : eps;
            x[i] = tmp+delta;       // bump this element
            delta = x[i] - tmp;     // try this to reduce error (from Schnabel)
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
        gelimd(jac,p,x,n);
    }
// Delete temporary storage
    delete [] p;
    for (i=0;i<n;i++) {
        delete [] jac[i];
    }
    delete [] jac;
}
/* This version executes up to 'maxiter' loops, with early exit
 * if the solution vector converges.
 */
void nlnewt3(void (*f)(double *x,double *fv,int n),double *x,
    double *fv,int n,double eps,int *maxiter)
{
    double tmp,delta,**jac,*p,*x0;
    int i,j,k;

    jac = new double *[n];
    for (i=0;i<n;i++) {
        jac[i] = new double [n];
    }
    p = new double [n];
    x0 = new double [n];

    for (k=0;k<*maxiter;k++) {
        f(x,fv,n);                  // get residuals for current value of 'x'

// Compute Jacobian matrix
        for (i=0;i<n;i++) {
            tmp = x[i];
            delta = (tmp > 1.0) ? eps*tmp : eps;
            x[i] = tmp+delta;       // bump this element
            delta = x[i] - tmp;     // try this to reduce error (from Schnabel)
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
        gelimd(jac,p,x0,n);
// Test for convergence
        tmp = 0.0;
        for (i=0;i<n;i++) {
            tmp += fabs(x[i]-x0[i]);
            x[i] = x0[i];
        }
        if (tmp < 1e-4) break;
    }
    *maxiter = k;
// Delete temporary storage
    delete [] x0;
    delete [] p;
    for (i=0;i<n;i++) {
        delete [] jac[i];
    }
    delete [] jac;
}
