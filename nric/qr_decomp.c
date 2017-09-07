

#include <math.h>

#include "nric.h"


  void
  qr_decomp(double** a, int n, double* c, double* d, int* sing)
  {
    *sing = 0;
    for (int k = 1; k < n; ++k)
      {
	double scale = 0.0;
	for (int i = k; i <= n; ++i)
	scale = dmax(scale, fabs(a[i][k]));
	if (scale == 0.0)
	  {
	    *sing = 1;
	    c[k] = d[k] = 0.0;
	  }
	else
	  {
	    for (int i = k; i <= n; ++i)
	      a[i][k] /= scale;
	    double sum = 0.0;
	    for (int i = k; i <= n; ++i)
	      sum += dsqr(a[i][k]);
	    double sigma = dsign(sqrt(sum), a[k][k]);
	    a[k][k] += sigma;
	    c[k] = sigma * a[k][k];
	    d[k] = -scale * sigma;
	    for (int j = k + 1; j <= n; ++j)
	      {
		double sum = 0.0;
		for (int i = k; i <= n; ++i)
		  sum += a[i][k] * a[i][j];
		double tau = sum / c[k];
		for (int i = k; i <= n; ++i)
		  a[i][j] -= tau * a[i][k];
	      }
	  }
      }
    d[n] = a[n][n];
    if (d[n] == 0.0)
      *sing = 1;
  }


  void
  qr_backsub(double **a, int n, double *c, double *d, double *b)
  {
    for (int j = 1; j < n; ++j)
      {
	double sum = 0.0;
	for (int i = j; i <= n; ++i)
	  sum += a[i][j] * b[i];
	double tau = sum / c[j];
	for (int i = j; i <= n; ++i)
	  b[i] -= tau * a[i][j];
      }
    r_backsub(a, n, d, b);
  }


  void
  r_backsub(double **a, int n, double *d, double *b)
  {
    b[n] /= d[n];
    for (int i = n - 1; i >= 1; --i)
      {
	double sum = 0.0;
	for (int j = i + 1; j <= n; ++j)
	  sum += a[i][j] * b[j];
	b[i] = (b[i] - sum) / d[i];
      }
  }


  void
  qr_update(double **r, double **qt, int n, double *u, double *v)
  {
    int k = n;
    for (; k >= 1; --k)
      if (u[k])
	break;
    if (k < 1)
      k = 1;
    for (int i = k - 1; i >= 1; --i)
      {
	jacobi_rotate(r, qt, n, i, u[i], -u[i+1]);
	if (u[i] == 0.0)
	  u[i] = fabs(u[i+1]);
	else if (fabs(u[i]) > fabs(u[i+1]))
	  u[i] = fabs(u[i]) * sqrt(1.0 + dsqr(u[i+1] / u[i]));
	else
	  u[i] = fabs(u[i+1]) * sqrt(1.0 + dsqr(u[i] / u[i+1]));
      }
    for (int j = 1; j <= n; ++j)
      r[1][j] += u[1] * v[j];
    for (int i = 1; i < k; ++i)
      jacobi_rotate(r, qt, n, i, r[i][i], -r[i+1][i]);
  }


  void
  jacobi_rotate(double** r, double** qt, int n, int i, double a, double b)
  {
    double c, s;
    if (a == 0.0)
      {
	c = 0.0;
	s = dsgn(b);
      }
    else
      {
	double fact = dpythag(a, b);
	c = a / fact;
	s = b / fact;
      }

    for (int j = i; j <= n; ++j)
      {
	double y = r[i][j];
	double w = r[i+1][j];
	r[i][j] = c * y - s * w;
	r[i+1][j] = s * y + c * w;
      }

    for (int j = 1; j <= n; ++j)
      {
	double y = qt[i][j];
	double w = qt[i+1][j];
	qt[i][j] = c * y - s * w;
	qt[i+1][j] = s * y + c * w;
      }
  }


