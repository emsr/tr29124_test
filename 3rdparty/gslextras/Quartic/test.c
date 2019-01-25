
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_heapsort.h>

/* sort by Re(z) then by Im(z) */
static int
cmp_cplx(const double *a, const double *b)
{
  double r = a[0] - b[0];

  if (r == 0.0)
    {
      double t = a[1] - b[1];
	    return t < 0.0 ? -1 : t > 0.0 ? 1 : 0;
    }
  else if (r < 0.0)
    return -1;
  else
    return 1;
}

int
main (void)
{
  const double eps = 100.0 * GSL_DBL_EPSILON;

  gsl_ieee_env_setup ();

  /* Quartic with complex roots */
  {
    gsl_complex z0, z1, z2, z3;
    int n =
      gsl_poly_complex_solve_quartic (0.0,0.0,0.0,-81.0, &z0, &z1, &z2,
        			     &z3);
    gsl_test (n != 4, "four roots, x^4-81");
    gsl_test_rel (GSL_REAL (z0), -3.0, 1e-9, "z0.real");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag");
    gsl_test_rel (GSL_REAL (z1), 0.0, 1e-9, "z1.real");
    gsl_test_rel (GSL_IMAG (z1), -3.0, 1e-9, "z1.imag");
    gsl_test_rel (GSL_REAL (z2), 0.0, 1e-9, "z2.real");
    gsl_test_rel (GSL_IMAG (z2), 3.0, 1e-9, "z2.imag");
    gsl_test_rel (GSL_REAL (z3), 3.0, 1e-9, "z3.real");
    gsl_test_rel (GSL_IMAG (z3), 0.0, 1e-9, "z3.imag");
  }
  {
    double sol=3.0/sqrt(2.0);
    gsl_complex z0, z1, z2, z3;
    int n =
      gsl_poly_complex_solve_quartic (0.0,0.0,0.0,81.0, &z0, &z1, &z2,
        			     &z3);
    gsl_test (n != 4, "four roots, x^4+81");
    gsl_test_rel (GSL_REAL (z0), -sol, 1e-9, "z0.real");
    gsl_test_rel (GSL_IMAG (z0), -sol, 1e-9, "z0.imag");
    gsl_test_rel (GSL_REAL (z1), -sol, 1e-9, "z1.real");
    gsl_test_rel (GSL_IMAG (z1), sol, 1e-9, "z1.imag");
    gsl_test_rel (GSL_REAL (z2), sol, 1e-9, "z2.real");
    gsl_test_rel (GSL_IMAG (z2), -sol, 1e-9, "z2.imag");
    gsl_test_rel (GSL_REAL (z3), sol, 1e-9, "z3.real");
    gsl_test_rel (GSL_IMAG (z3), sol, 1e-9, "z3.imag");
  }
  {
    gsl_complex z0, z1, z2, z3;
    int n =
      gsl_poly_complex_solve_quartic (0.0,4.0,0.0,0.0, &z0, &z1, &z2,
        			     &z3);
    gsl_test (n != 4, "four roots, x^2*(x^2+4)");
    gsl_test_rel (GSL_REAL (z0), 0.0, 1e-9, "z0.real");
    gsl_test_rel (GSL_IMAG (z0), -2.0, 1e-9, "z0.imag");
    gsl_test_rel (GSL_REAL (z1), 0.0, 1e-9, "z1.real");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag");
    gsl_test_rel (GSL_REAL (z2), 0.0, 1e-9, "z2.real");
    gsl_test_rel (GSL_IMAG (z2), 0.0, 1e-9, "z2.imag");
    gsl_test_rel (GSL_REAL (z3), 0.0, 1e-9, "z3.real");
    gsl_test_rel (GSL_IMAG (z3), 2.0, 1e-9, "z3.imag");
  }
  {
    gsl_complex z0, z1, z2, z3;
    int n =
      gsl_poly_complex_solve_quartic (0.0,-4.0,0.0,0.0, &z0, &z1, &z2,
        			     &z3);
    gsl_test (n != 4, "four roots, x^2*(x-2)*(x+2)");
    gsl_test_rel (GSL_REAL (z0), -2.0, 1e-9, "z0.real");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag");
    gsl_test_rel (GSL_REAL (z1), 0.0, 1e-9, "z1.real");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag");
    gsl_test_rel (GSL_REAL (z2), 0.0, 1e-9, "z2.real");
    gsl_test_rel (GSL_IMAG (z2), 0.0, 1e-9, "z2.imag");
    gsl_test_rel (GSL_REAL (z3), 2.0, 1e-9, "z3.real");
    gsl_test_rel (GSL_IMAG (z3), 0.0, 1e-9, "z3.imag");
  }
  {
    gsl_complex z0, z1, z2, z3;
    int n =
      gsl_poly_complex_solve_quartic (-10.0,35.0,-50.0,24.0, &z0, &z1, &z2,
        			     &z3);
    gsl_test (n != 4, "four roots, (x-1)*(x-2)*(x-3)*(x-4)");
    gsl_test_rel (GSL_REAL (z0), 1.0, 1e-9, "z0.real");
    gsl_test_rel (GSL_IMAG (z0), 0.0, 1e-9, "z0.imag");
    gsl_test_rel (GSL_REAL (z1), 2.0, 1e-9, "z1.real");
    gsl_test_rel (GSL_IMAG (z1), 0.0, 1e-9, "z1.imag");
    gsl_test_rel (GSL_REAL (z2), 3.0, 1e-9, "z2.real");
    gsl_test_rel (GSL_IMAG (z2), 0.0, 1e-9, "z2.imag");
    gsl_test_rel (GSL_REAL (z3), 4.0, 1e-9, "z3.real");
    gsl_test_rel (GSL_IMAG (z3), 0.0, 1e-9, "z3.imag");
  }
  {
    gsl_complex z0, z1, z2, z3;
    int n =
      gsl_poly_complex_solve_quartic (-5.0,7.0,-5.0,6.0, &z0, &z1, &z2,
        			     &z3);
    gsl_test (n != 4, "four roots, (x-i)*(x+i)*(x-2)*(x-3)");
    gsl_test_rel (GSL_REAL (z0), 0.0, 1e-9, "z0.real");
    gsl_test_rel (GSL_IMAG (z0), -1.0, 1e-9, "z0.imag");
    gsl_test_rel (GSL_REAL (z1), 0.0, 1e-9, "z1.real");
    gsl_test_rel (GSL_IMAG (z1), 1.0, 1e-9, "z1.imag");
    gsl_test_rel (GSL_REAL (z2), 2.0, 1e-9, "z2.real");
    gsl_test_rel (GSL_IMAG (z2), 0.0, 1e-9, "z2.imag");
    gsl_test_rel (GSL_REAL (z3), 3.0, 1e-9, "z3.real");
    gsl_test_rel (GSL_IMAG (z3), 0.0, 1e-9, "z3.imag");
  }
  {
    gsl_complex z0, z1, z2, z3;
    int n =
      gsl_poly_complex_solve_quartic (0.0, 5.0, 0.0, 4.0, &z0, &z1, &z2, &z3);
    gsl_test (n != 4, "four roots, (x-i)(x+i)(x-2i)(x+2i)=0");
    gsl_test_rel (GSL_REAL (z0), 0.0, 1e-9,"z0.real");
    gsl_test_rel (GSL_IMAG (z0), -2.0, 1e-9,"z0.imag");
    gsl_test_rel (GSL_REAL (z1), 0.0, 1e-9,"z1.real");
    gsl_test_rel (GSL_IMAG (z1), -1.0, 1e-9,"z1.imag");
    gsl_test_rel (GSL_REAL (z2), 0.0, 1e-9,"z2.real");
    gsl_test_rel (GSL_IMAG (z2), 1.0, 1e-9,"z2.imag");
    gsl_test_rel (GSL_REAL (z3), 0.0, 1e-9,"z3.real");
    gsl_test_rel (GSL_IMAG (z3), 2.0, 1e-9,"z3.imag");
  }
  {
    gsl_complex z0, z1, z2, z3;
    int n =
      gsl_poly_complex_solve_quartic (-4.0, 11.0, -14.0, 10.0, &z0, &z1, &z2,
        			     &z3);
    gsl_test (n != 4, "four roots, (x-1-i)(x-1+i)(x-1-2i)(x-1+2i)=0");
    gsl_test_rel (GSL_REAL (z0), 1.0, 1e-9,"z0.real");
    gsl_test_rel (GSL_IMAG (z0), -2.0, 1e-9,"z0.imag");
    gsl_test_rel (GSL_REAL (z1), 1.0, 1e-9,"z1.real");
    gsl_test_rel (GSL_IMAG (z1), -1.0, 1e-9,"z1.imag");
    gsl_test_rel (GSL_REAL (z2), 1.0, 1e-9,"z2.real");
    gsl_test_rel (GSL_IMAG (z2), 1.0, 1e-9,"z2.imag");
    gsl_test_rel (GSL_REAL (z3), 1.0, 1e-9,"z3.real");
    gsl_test_rel (GSL_IMAG (z3), 2.0, 1e-9,"z3.imag");
  }


  {
    double x0, x1, x2, x3;
    int n =
      gsl_poly_solve_quartic (0.0,0.0,0.0,-81.0,&x0,&x1,&x2,&x3);
    gsl_test (n != 2, "two real roots, x^4-81");
    gsl_test_rel (x0, -3.0, 1e-9, "x0");
    gsl_test_rel (x1, 3.0, 1e-9, "x1");
  }
  {
    double x0, x1, x2, x3;
    int n =
      gsl_poly_solve_quartic (0.0,0.0,0.0,81.0,&x0,&x1,&x2,&x3);
    gsl_test (n != 0, "zero roots, x^4+81");
  }
  {
    double x0, x1, x2, x3;
    int n =
      gsl_poly_solve_quartic (0.0,-4.0,0.0,0.0,&x0,&x1,&x2,&x3);
    gsl_test (n != 4, "four roots, x^2*(x-2)*(x+2)");
    gsl_test_rel (x0, -2.0, 1e-9, "x0");
    gsl_test_rel (x1, 0.0, 1e-9, "x1");
    gsl_test_rel (x2, 0.0, 1e-9, "x2");
    gsl_test_rel (x3, 2.0, 1e-9, "x3");
  }
  {
    double x0, x1, x2, x3;
    int n =
      gsl_poly_solve_quartic (-10.0,35.0,-50.0,24.0,&x0,&x1,&x2,&x3);
    gsl_test (n != 4, "four roots, (x-1)*(x-2)*(x-3)*(x-4)");
    gsl_test_rel (x0, 1.0, 1e-9, "x0");
    gsl_test_rel (x1, 2.0, 1e-9, "x1");
    gsl_test_rel (x2, 3.0, 1e-9, "x2");
    gsl_test_rel (x3, 4.0, 1e-9, "x3");
  }

  /* now summarize the results */

  exit (gsl_test_summary ());
}
