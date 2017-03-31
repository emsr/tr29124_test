/*  cged.c -- Gaussian elimination linear equation solver for complex
 *            double precision numbers.
 *  (C)2001, C. Bond. All rights reserved.
 *
 *      Simple pivoting on zero diagonal entry supported.
 *      This version eliminates unnecessary zeroing of lower
 *      triangular portion of source matrix. Does not scale
 *      rows to unity pivot value. Swaps b[] as well as a[][],
 *      so pivot ID vector is not required.
 */
#include <math.h>
#include <complex.h>

int cged(complex<double> **a,complex<double> *b,
    complex<double> *x, int n)
{
    complex<double> tmp,pvt,*t;
    int i,j,k,itmp;

    for (i=0;i<n;i++) {             // outer loop on rows
        pvt = a[i][i];              // get pivot value
        if (!abs(pvt)) {
            for (j=i+1;j<n;j++) {
                if ((pvt = abs(a[j][i])) != 0.0) break;
            }
            if (!abs(pvt))
                return 1;           // nowhere to run!
            t=a[j];                 // swap matrix rows...
            a[j]=a[i];
            a[i]=t;
            tmp=b[j];               // ...and result vector
            b[j]=b[i];
            b[i]=tmp;        
        }
        for (k=i+1;k<n;k++) {       // (virtual) eliminate column
            tmp = a[k][i]/pvt;
            for (j=i+1;j<n;j++) {
                a[k][j] -= tmp*a[i][j];
            }
            b[k] -= tmp*b[i];
		}
	}
// Do back substitution
	for (i=n-1;i>=0;i--) {
		x[i]=b[i];
		for (j=n-1;j>i;j--) {
			x[i] -= a[i][j]*x[j];
		}
		x[i] /= a[i][i];
	}
    return 0;
}
/* This version preserves the matrix 'a' and the vector 'b'. */

int cged2(complex<double> **a,complex<double> *b,
    complex<double> *x, int n)
{
    complex<double> tmp,pvt,*t,**aa,*bb;
    int i,j,k,itmp,retval;

// Initialize return value for successful execution.
    retval = 0;

// Create and initialize working storage
    aa = new complex<double> *[n];
    bb = new complex<double> [n];
    for (i=0;i<n;i++) {
        aa[i] = new complex<double> [n];
        for (j=0;j<n;j++) {
            aa[i][j] = a[i][j];
        }
        bb[i] = b[i];
    }

// Main loop
    for (i=0;i<n;i++) {             // outer loop on rows
        pvt = aa[i][i];              // get pivot value
        if (!abs(pvt)) {
            for (j=i+1;j<n;j++) {
                if ((pvt = abs(aa[j][i])) != 0.0) break;
            }
            if (!abs(pvt)) {
                retval = 1;
                goto _100;          // pull the plug!
            }
            t=aa[j];                 // swap matrix rows...
            aa[j]=aa[i];
            aa[i]=t;
            tmp=bb[j];               // ...and result vector
            bb[j]=bb[i];
            bb[i]=tmp;        
        }
        for (k=i+1;k<n;k++) {       // (virtual) eliminate column
            tmp = aa[k][i]/pvt;
            for (j=i+1;j<n;j++) {
                aa[k][j] -= tmp*aa[i][j];
            }
            bb[k] -= tmp*bb[i];
		}
	}
// Do back substitution
	for (i=n-1;i>=0;i--) {
        x[i]=bb[i];
		for (j=n-1;j>i;j--) {
            x[i] -= aa[i][j]*x[j];
		}
        x[i] /= aa[i][i];
	}
_100:
    for (i=0;i<n;i++) {
        delete [] aa[i];
    }
    delete [] aa;
    delete [] bb;
    return retval;
}
