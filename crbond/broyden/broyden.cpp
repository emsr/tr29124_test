/*  broyden.cpp -- Broyden's method for solving a nonlinear system of
 *              equations.
 *  (C) 2001, C. Bond. All rights reserved.
 */
#include <math.h>

int gelimd2(double **a,double *b,double *x,int n);


int broyden(void (*f)(double *x,double *fv,int n),
    double *x0,double *f0,int n,double *eps,int *iter)
{
    double **A,*x1,*f1,*s,d,tmp;
    int i,j,k;

// Allocate temporary memory
    x1 = new double [n];
    f1 = new double [n];
    s = new double [n];
    A = new double *[n];
    for (i=0;i<n;i++) {
        A[i] = new double [n];
    }

// Create identity matrix for startup (consider Jacobian alternative)
    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) {
            A[i][j] = 0.0;
        }
        A[i][i] = 1.0;
    }
// Main loop
    for (k=0;k< *iter;k++){
        f(x0,f0,n);
        gelimd2(A,f0,s,n);      // Signs of 'f0', 's' are reversed
        d = 0.0;
        for (i=0;i<n;i++) {
            x1[i] = x0[i] - s[i];
            d += s[i]*s[i];
        }
        if (d <= fabs(*eps)) break;
        f(x1,f1,n);
// Update A
        for (i=0;i<n;i++) {
            tmp = f1[i]+f0[i];
            for (j=0;j<n;j++) {
                tmp -= A[i][j]*s[j];
            }
            for (j=0;j<n;j++) {
                A[i][j] -= tmp*s[j]/d;
            }
        }
        for (i=0;i<n;i++) {
            f0[i] = f1[i];
            x0[i] = x1[i];
        }
    }
    *iter = k;
// Deallocate memory
    for (i=0;i<n;i++) {
        delete [] A[i];
    }
    delete [] A;
    delete [] x1;
    delete [] f1;
    delete [] s;
    return 0;
}
