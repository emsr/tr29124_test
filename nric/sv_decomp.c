

#include <math.h>


#include "nric.h"


  void
  sv_decomp(double** a, int m, int n, double *w, double **v)
  {
    int flag, i, its, j, jj, k, l, nm;
    double c, f, h, s, x, y, z;

    const int ITS = 30;

    double* rv1 = dvector(1,n);

    double g = 0.0;
    double scale = 0.0;
    double anorm = 0.0;
    /*
     * Householder reduction to bidiagonal form.
     */
    for (int i = 1; i <= n; ++i)
      {
	l = i + 1;
	rv1[i] = scale * g;
	g = s = scale = 0.0;
	if (i <= m)
	  {
	    for (int k = i; k <= m; ++k)
	      scale += fabs(a[k][i]);
	    if (scale)
	      {
		for (int k = i; k <= m; ++k)
		  {
		    a[k][i] /= scale;
		    s += a[k][i] * a[k][i];
		  }
		double f = a[i][i];
		g = -dsign(sqrt(s), f);
		double h = f * g - s;
		a[i][i] = f - g;
		for (int j = l; j <= n; ++j)
		  {
		    for (s = 0.0, k = i; k <= m; ++k)
		      s += a[k][i] * a[k][j];
		    f = s / h;
		    for (k = i; k <= m; ++k)
		      a[k][j] += f * a[k][i];
		  }
		for (int k = i; k <= m; ++k)
		  a[k][i] *= scale;
	      }
	  }
	w[i] = scale * g;
	g = s = scale = 0.0;
	if (i <= m && i != n)
	  {
	    for (int k = l; k <= n; ++k)
	      scale += fabs(a[i][k]);
	    if (scale)
	      {
		for (int k = l; k <= n; ++k)
		  {
		    a[i][k] /= scale;
		    s += a[i][k] * a[i][k];
		  }
		double f = a[i][l];
		g = -dsign(sqrt(s), f);
		double h = f * g - s;
		a[i][l] = f - g;
		for (int k = l; k <= n; ++k)
		  rv1[k] = a[i][k]/h;
		for (int j = l; j <= m; ++j)
		  {
		    double s = 0.0;
		    for (int k = l; k <= n; ++k)
		      s += a[j][k] * a[i][k];
		    for (int k = l; k <= n; ++k)
		      a[j][k] += s * rv1[k];
		  }
		for (k = l; k <= n; ++k)
		  a[i][k] *= scale;
	      }
	  }
	anorm = dmax(anorm, fabs(w[i]) + fabs(rv1[i]));
      }

    /*
     * Accumulation of right-hand decomposition V.
     */
    for (int i = n; i >= 1; --i)
      {
	if (i < n)
	  {
	    if (g)
	      {
		for (int j = l; j <= n; ++j)
		  v[j][i] = (a[i][j] / a[i][l]) / g;
		for (int j = l; j <= n; ++j)
		  {
		    double s = 0.0;
		    for (int k = l; k <= n; ++k)
		      s += a[i][k] * v[k][j];
		    for (int k = l; k <= n; ++k)
		      v[k][j] += s * v[k][i];
		  }
	      }
	    for (j = l; j <= n; ++j)
	      v[i][j] = v[j][i] = 0.0;
	  }
	v[i][i] = 1.0;
	g = rv1[i];
	l = i;
      }

    /*
     * Accumulation of left-hand decompositions.
     */
    for (i = imin(m, n); i >= 1; --i)
      {
	l = i + 1;
	g = w[i];
	for (int j = l; j <= n; ++j)
	  a[i][j] = 0.0;
	if (g)
	  {
	    g = 1.0/g;
	    for (int j = l; j <= n; ++j)
	      {
		double s = 0.0;
		for (int k = l; k <= m; ++k)
		  s += a[k][i] * a[k][j];
		f = (s / a[i][i]) * g;
		for (int k = i; k <= m; ++k)
		  a[k][j] += f * a[k][i];
	      }
	    for (int j = i; j <= m; ++j)
	      a[j][i] *= g;
	  }
	else
	  for (int j = i; j <= m; ++j)
	    a[j][i] = 0.0;
	++a[i][i];
      }

    /*
     * Diagonalization of the bidiagonal form;
     */
    for (k = n; k >= 1; --k)
      {
	for (its = 1; its <= ITS; ++its)
	  {
	    flag = 1;
	    for (l = k; l >= 1; --l)
	      {
		nm = l - 1;
		if (fabs(rv1[l]) + anorm == anorm)
		  {
		    flag = 0;
		    break;
		  }
		if (fabs(w[nm]) + anorm == anorm)
		  break;
	      }
	    if (flag)
	      {
		c = 0.0;
		s = 1.0;
		for (i = l; i <= k; ++i)
		  {
		    f = s * rv1[i];
		    rv1[i] = c * rv1[i];
		    if (fabs(f) + anorm == anorm)
		      break;
		    g = w[i];
		    h = dpythag(f, g);
		    w[i] = h;
		    h = 1.0 / h;
		    c = g * h;
		    s = -f * h;
		    for (j = 1; j <= m; ++j)
		      {
			y = a[j][nm];
			z = a[j][i];
			a[j][nm] = y * c + z * s;
			a[j][i] = z * c - y * s;
		      }
		  }
	      }
	    z = w[k];
	    if (l == k)
	      {
		/*
		 * Convergence!!!
		 */
		if (z < 0.0)
		  {
		    /*
		     * Make singular value non negative.
		     */
		    w[k] = -z;
		    for (j = 1; j <= n; ++j)
		      v[j][k] = -v[j][k];
		  }
		break;
	      }
	    if (its == ITS)
		nrerror("No convergence in  30 sv_decomp iterations.");

	    /*
	     * Shift from bottom 2x2 minor.
	     */
	    x = w[l];
	    nm = k - 1;
	    y = w[nm];
	    g = rv1[nm];
	    h = rv1[k];
	    f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
	    g = dpythag(f, 1.0);
	    f = ((x - z) * (x + z) + h * ((y / (f + dsign(g, f))) - h)) / x;

	    /*
	     * Next QR transformation.
	     */
	    c = s = 1.0;
	    for (j = l; j <= nm; ++j)
	      {
		i = j + 1;
		g = rv1[i];
		y = w[i];
		h = s * g;
		g = c * g;
		z = dpythag(f, h);
		rv1[j] = z;
		c = f / z;
		s = h / z;
		f = x * c + g * s;
		g = g * c - x * s;
		h = y * s;
		y *= c;
		for (jj = 1; jj <= n; ++jj)
		  {
		    x = v[jj][j];
		    z = v[jj][i];
		    v[jj][j] = x * c + z * s;
		    v[jj][i] = z * c - x * s;
		  }
		z = dpythag(f,h);
		w[j] = z;
		/*
		 * Rotation can be arbitrary if z = 0.
		 */
		if (z)
		  {
		    z = 1.0 / z;
		    c = f * z;
		    s = h * z;
		  }
		f = c * g + s * y;
		x = c * y - s * g;
		for (jj = 1; jj <= m; ++jj)
		  {
		    y = a[jj][j];
		    z = a[jj][i];
		    a[jj][j] = y*c + z*s;
		    a[jj][i] = z*c - y*s;
		  }
	      }
	    rv1[l] = 0.0;
	    rv1[k] = f;
	    w[k] = x;
	  }
      }
    free_dvector(rv1, 1, n);
  }




  void
  sv_backsub(double** u, double* w, double** v, int m, int n, double *b, double *x)
  {
    double* tmp = dvector(1, n);
    for (int j = 1; j <= n; ++j)
      {
	double s = 0.0;
	if (w[j])
	  {
	    for (int i = 1; i <= m; ++i)
	      s += u[i][j] * b[i];
	    s /= w[j];
	  }
	tmp[j] = s;
      }
    for (int j = 1; j <=n; ++j)
      {
	double s = 0.0;
	for (int jj = 1; jj <=n; ++jj)
	  s += v[j][jj] * tmp[jj];
	x[j] = s;
      }
    free_dvector(tmp, 1, n);
  }



  /*
   *    Improves a solution vector x of the linear set A.x = b.  The matrix a and the
   *    SV decomposition of a -- u, w, v and the
   *    right-hand side vector are input along with the solution vector x.
   *    The solution vector x is improved and modified on output.
   */
  void
  sv_improve(double** a, double** u, double* w, double** v, int m, int n, double* b, double* x)
  {
    double* r = dvector(1, m);
    double* dx = dvector(1, n);

    for (int i = 1; i <= m; ++i)
      {
	r[i] = -b[i];
	for (int j = 1; j <= n; ++j)
	  r[i] += a[i][j]*x[j];
      }
    sv_backsub(u, w, v, m, n, r, dx);
    for (int i = 1; i <= n; ++i)
      x[i] -= dx[i];

    free_dvector(dx, 1, n);
    free_dvector(r, 1, m);
  }



