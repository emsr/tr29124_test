/* ged.cpp -- Gaussian elimination linear equation solvers.
 *
 *  (C) 2001, C. Bond. All rights reserved.
 *
 *  Simple pivoting on zero diagonal element supported.
 *   Eliminates unnecessary  zeroing of lower triangle.
 *   Does not scale rows to unity pivot value.
 *   Swaps b[] as well as a[][], so a pivot ID vector
 *   is not required.
 */
#include <math.h>

int gelimd(double **a,double *b,double *x, int n)
{
    double tmp,pvt,*t;
    int i,j,k,itmp;

    for (i=0;i<n;i++) {             // outer loop on rows
        pvt = a[i][i];              // get pivot value
        if (!pvt) {
            for (j=i+1;j<n;j++) {
                if((pvt = a[j][i]) != 0.0) break;
            }
            if (!pvt) return 1;     // nowhere to run!
            t=a[j];                 // swap matrix rows...
            a[j]=a[i];
            a[i]=t;
            tmp=b[j];               // ...and result vector
            b[j]=b[i];
            b[i]=tmp;        
        }
// (virtual) Gaussian elimination of column
        for (k=i+1;k<n;k++) {       // alt: for (k=n-1;k>i;k--)
            tmp = a[k][i]/pvt;
            for (j=i+1;j<n;j++) {   // alt: for (j=n-1;j>i;j--)
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
int gelimd2(double **a,double *b,double *x, int n)
{
    double tmp,pvt,*t,**aa,*bb;
    int i,j,k,itmp,retval;

    retval = 0;
    aa = new double *[n];
    bb = new double [n];
    for (i=0;i<n;i++) {
        aa[i] = new double [n];
        bb[i] = b[i];
        for (j=0;j<n;j++) {
            aa[i][j] = a[i][j];
        }
    }
    for (i=0;i<n;i++) {             // outer loop on rows
        pvt = aa[i][i];              // get pivot value
        if (!pvt) {
            for (j=i+1;j<n;j++) {
                if((pvt = a[j][i]) != 0.0) break;
            }
            if (!pvt) {
                retval = 1;
                goto _100;     // pull the plug!
            }
            t=aa[j];                 // swap matrix rows...
            aa[j]=aa[i];
            aa[i]=t;
            tmp=bb[j];               // ...and result vector
            bb[j]=bb[i];
            bb[i]=tmp;        
        }
// (virtual) Gaussian elimination of column
        for (k=i+1;k<n;k++) {       // alt: for (k=n-1;k>i;k--)
            tmp = aa[k][i]/pvt;
            for (j=i+1;j<n;j++) {   // alt: for (j=n-1;j>i;j--)
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
// Deallocate memory
_100:
    for (i=0;i<n;i++) {
        delete [] aa[i];
    }
    delete [] aa;
    delete [] bb;
    return retval;
}

