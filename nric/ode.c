

#include <math.h>

#include "nric.h"




/*
 *    Given values for n dependent variables y[1..n] and their derivatives dydx[1..n] known at x,
 *    use fourth-order Runge-Kutta method to advance the solution over an interval h and returns
 *    the incremented variables as yout[1..n], (which need not be a distinct array from y).
 *    The user supplies the routine derivs which returns the derivatives dydx[1..n] at x.
 */
void
runge_kutta_4(double *y, double *dydx, int n,
              double x, double h, double *yout,
              void (*derivs)(double, double *, double *))
{
  double * dym = dvector(1, n);
  double * dyt = dvector(1, n);
  double * yt = dvector(1, n);

  double h2 = h / 2.0;
  double h6 = h / 6.0;
  double xh = x + h2;
  for (int i = 1; i <= n; ++i)
    yt[i] = y[i] + h2*dydx[i];
  (*derivs)(xh, yt, dyt);
  for (int i = 1; i <= n; ++i)
    yt[i] = y[i] + h2*dyt[i];
  (*derivs)(xh, yt, dym);
  for (int i = 1; i <= n; ++i)
    {
      yt[i] = y[i] + h*dym[i];
      dym[i] += dyt[i];
    }
  (*derivs)(x+h, yt, dyt);
  for (int i = 1; i <= n; ++i)
    yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);

  free_dvector(yt, 1, n);
  free_dvector(dyt, 1, n);
  free_dvector(dym, 1, n);
}


/*
 *    Variables for the storage of Runge - Kutta intermediate results.
 *    You have to allocate space for these things:
 *    runge_kutta_xx[1..nvar] and runge_kutta_yy[1..nstep+1][1..nstep+1]
 *    where nvar is the number of dependent variables and nstep is the number of steps.
 *    These arrays are set in dumb_runge_kutta.
 */


/*
 *    Starting from initial values ystart[1..nvar] for nvar functions known at x1,
 *    use 4th order Runge-Kutta to advance nstep equal increments to x2.
 *    The user supplied routine derivs supplies derivatives.
 *    The results are stored in the global variables runge_kutta_xx[1..nstep+1], runge_kutta_yy[1..nvar][1..nstep+1].
 */
void
dumb_runge_kutta(double *ystart, int nvar, double x1, double x2, int nstep,
                 double *runge_kutta_xx, double **runge_kutta_yy,
                 void (*derivs)(double, double *, double *))
{
  double * y = dvector(1, nvar);
  double * yout = dvector(1, nvar);
  double * dy = dvector(1, nvar);

  //  Load starting values of dependant variables.
  for (int i = 1; i <= nvar; ++i)
    {
      y[i] = ystart[i];
      runge_kutta_yy[i][1] = y[i];
    }
  runge_kutta_xx[1] = x1;
  double x = x1;
  double h = (x2 - x1) / nstep;

  //  Take nstep steps in the independant variable x.
  for (int k = 1; k <= nstep; ++k)
    {
      (*derivs)(x, y, dy);
      runge_kutta_4(y, dy, nvar, x, h, yout, derivs);
      if ((x + h) == x)
        nrerror("Step size too small in dumb_runge_kutta.");
      x += h;
      //  Store intermediate results.
      runge_kutta_xx[k+1] = x;
      for (int i = 1; i <= nvar; ++i)
        {
          y[i] = yout[i];
          runge_kutta_yy[i][k + 1] = y[i];
        }
    }

  free_dvector(dy, 1, nvar);
  free_dvector(yout, 1, nvar);
  free_dvector(y, 1, nvar);
}




/*
 *    Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure
 *    accuracy and adjust stepsize.  Input are the dependent variable vector y[1..n]
 *    and its derivative dydx[1..n] at the starting value of the independent variable x.
 *    Also input are the first guess for the stepsize htry, the requred accuracy eps, and the
 *    vector yscale[1..n] against which the error is scaled independently for each dependent variable.
 *    On output, x and y are replaced by thier new values, hdid is the stepsize which was actually
 *    accomplished, and hnext is the estimated next stepsize.
 *    derivs is the user-supplied routine which computes the right-hand side derivatives.
 */

void
quad_runge_kutta(double *y, double *dydx, int n, double *x, double htry,
                 double eps, double *yscale, double *hdid, double *hnext,
                 void (*derivs)(double, double *, double *))
{
  const double POW_GROW = 10000;
  const double POW_SHRINK = 1.0e-30;
  const double F_CORR = 1.0/15.0;
  const double F_SAFETY = 0.9;
  const double ERR_COND = pow((4.0/F_SAFETY), (1.0/POW_GROW));

  double * dysav = dvector(1, n);
  double * ysav = dvector(1, n);
  double * ytemp = dvector(1, n);
  double xsav = *x;
  for (int i = 1; i <= n; ++i)
    {

      ysav[i] = y[i];
      dysav[i] = dydx[i];
    }
  //  Set stepsize to the initial trial value.
  double h = htry;
  while (1)
    {

      //  Take two half steps.
      double h2 = h / 2.0;
      runge_kutta_4(ysav, dysav, n, xsav, h2, ytemp, derivs);
      *x = xsav + h2;
      (*derivs)(*x, ytemp, dydx);
      runge_kutta_4(ytemp, dydx, n, *x, h2, y, derivs);
      *x = xsav + h;
      if (*x == xsav)
        nrerror("Step size too small in quad_runge_kutta.");
      //  Take the large step.
      runge_kutta_4(ysav, dysav, n, xsav, h, ytemp, derivs);
      //  Evaluate accuracy.  Put the error estimate into ytemp.
      double errmax = 0.0;
      for (int i = 1; i <= n; ++i)
        {
          ytemp[i] = y[i] - ytemp[i];
          double temp = fabs(ytemp[i]/yscale[i]);
          if (errmax < temp)
            errmax = temp;
        }
      //  Scale relative to required tolerance.
      errmax /= eps;
      if (errmax <= 1.0)
        {
          //  Step succeeded.  Compute size of next step.
          *hdid = h;
          *hnext = (errmax > ERR_COND ? F_SAFETY * exp(POW_GROW * log(errmax)) : 4.0 * h);
          break;
        }
      //  Truncation error too large, reduce stepsize.
      h = F_SAFETY * h * exp(POW_SHRINK * log(errmax));
    }
  //  Mop up fifth order truncation error.
  for (int i = 1; i <= n; ++i)
    y[i] += F_CORR * ytemp[i];

  free_dvector(ytemp, 1, n);
  free_dvector(ysav, 1, n);
  free_dvector(dysav, 1, n);
  return;
}



/*
 *
 */
void
cash_karp_rk(double *y, double *dydx, int n,
             double x, double h, double *yout, double *yerr,
             void (*derivs)(double, double *, double *))
{
    static double
      a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
      b21 = 0.2,
      b31 = 3.0/40.0, b32 = 9.0/40.0,
      b41 = 0.3, b42 = -0.9, b43 = 1.2,
      b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0, b54 = 35.0/27.0,
      b61 = 1631.0/55296.0, b62 = 175.0/512.0, b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
      c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
      dc5 = -277.0/14336.0;
    double dc1 = c1 - 2825.0/27648.0, dc3 = c3 - 18575.0/48384.0,
           dc4 = c4 - 13525.0/55296.0, dc6 = c6 - 0.25;

    double * ak2 = dvector(1, n);
    double * ak3 = dvector(1, n);
    double * ak4 = dvector(1, n);
    double * ak5 = dvector(1, n);
    double * ak6 = dvector(1, n);
    double * ytemp = dvector(1, n);

    for (int i = 1; i <= n; ++i)
      ytemp[i] = y[i] + h * b21 * dydx[i];
    (*derivs)(x + a2 * h, ytemp, ak2);
    for (int i = 1; i <= n; ++i)
      ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
    (*derivs)(x + a3 * h, ytemp, ak3);
    for (int i = 1; i <= n; ++i)
      ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
    (*derivs)(x + a4 * h, ytemp, ak4);
    for (int i = 1; i <= n; ++i)
      ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
    (*derivs)(x + a5 * h, ytemp, ak5);
    for (int i = 1; i <= n; ++i)
      ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
    (*derivs)(x + a6 * h, ytemp, ak6);
    for (int i = 1; i <= n; ++i)
      yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
    for (int i = 1; i <= n; ++i)
      yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);

    free_dvector(ytemp, 1, n);
    free_dvector(ak6, 1, n);
    free_dvector(ak5, 1, n);
    free_dvector(ak4, 1, n);
    free_dvector(ak3, 1, n);
    free_dvector(ak2, 1, n);
}




/*
 *    Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure
 *    accuracy and adjust stepsize.  Input are the dependent variable vector y[1..n]
 *    and its derivative dydx[1..n] at the starting value of the independent variable x.
 *    Also input are the first guess for the stepsize htry, the requred accuracy eps, and the
 *    vector yscale[1..n] against which the error is scaled independently for each dependent variable.
 *    On output, x and y are replaced by thier new values, hdid is the stepsize which was actually
 *    accomplished, and hnext is the estimated next stepsize.
 *    derivs is the user-supplied routine which computes the right-hand side derivatives.
 */

void
quad_cash_karp_rk(double *y, double *dydx, int n, double *x, double htry,
                  double eps, double *yscale, double *hdid, double *hnext,
                  void (*derivs)(double, double *, double *))
{
  int i;
  double errmax, h, htemp, xnew, *yerr, *ytemp;

  const double POW_GROW = -0.20;
  const double POW_SHRINK = -0.25;
  const double F_SAFETY = 0.9;
  const double ERR_COND = pow(5.0/F_SAFETY, 1.0/POW_GROW);

  yerr = dvector(1, n);
  ytemp = dvector(1, n);

  h = htry;
  while (1)
    {
      cash_karp_rk(y, dydx, n, *x, h, ytemp, yerr, derivs);
      errmax = 0.0;
      for (int i = 1; i <= n; ++i)
        errmax = dmax(errmax, fabs(yerr[i]/yscale[i]));
      errmax /= eps;
      if (errmax <= 1.0)
        break;
      htemp = F_SAFETY * h * pow(errmax, POW_SHRINK);
      h = (h > 0.0 ? dmax(htemp, 0.1 * h) : dmin(htemp, 0.1 * h));
      xnew= *x + h;
      if (xnew == *x)
        nrerror("Stepsize underflow in quad_cash_karp_rk.");
    }

  if (errmax > ERR_COND)
    *hnext = F_SAFETY * h * pow(errmax, POW_GROW);
  else
    *hnext = 5.0 * h;
  *x += *hdid = h;
  for (int i = 1; i <= n; ++i)
    y[i] = ytemp[i];

  free_dvector(ytemp, 1, n);
  free_dvector(yerr, 1, n);
}




/*
 *    Global variables for the storage of intermediate results from ode_integrate.
 */
int ode_max = 0;
int ode_count = 0;
double *ode_xp = 0;
double **ode_yp = 0;
double ode_dxsave = 0;


/*
 *    ODE driver with adaptive stepsize control.  Integrate starting with values
 *    ystart[1..nvar] from x1 to x2 with accuracy eps, storing intermediate results
 *    in global variables ode_xp, ode_yp, ode_max, ode_count, ode_dxsave.  If ode_max == 0 no intermediate results
 *    will be stored and the pointers ode_xp and ode_yp can be set to zero.
 *    h1 should be set as a first guess initial stepsize, hmin is the minimum stepsize (can be zero).
 *    On output nok and nbad are the numbers of good and bad (but retried and fixed) steps taken.
 *    ystart is replaced by stepped values at the end of the integration interval.
 *    derivs is a user-supplied function for returning the right side derivatives of the independent variables.
 *    stepper is the name of the integration stepper to be used (e.g. quad_runge_kutta or bulirsch_stoer).
 */
void
ode_integrate(double *ystart, int nvar, double x1, double x2,
              double eps, double h1, double hmin,
              int *nok, int *nbad,
              void (*derivs)(double, double *, double *),
              void (*stepper)(double *, double *, int, double *, double,
                              double, double *, double *, double *,
                              void (*)(double, double *, double *)))
{
  double xsave, hnext, hdid;

  const int MAXSTEP = 10000;
  const double TINY = 1.0e-30;

  double *yscale = dvector(1, nvar);
  double *y = dvector(1, nvar);
  double *dydx = dvector(1, nvar);

  double x = x1;
  double h = (x2 > x1) ? fabs(h1) : -fabs(h1);
  *nok = *nbad = ode_count = 0;
  for (int i = 1; i <= nvar; ++i)
    y[i] = ystart[i];
  if (ode_max > 0)
    xsave = x - 2.0 * ode_dxsave;
  for (int nstep = 1; nstep <= MAXSTEP; ++nstep)
   {
      (*derivs)(x, y, dydx);
      for (int i = 1; i <= nvar; ++i)
        yscale[i] = fabs(y[i]) + fabs(dydx[i]) + TINY;
      if (ode_max)
        if (fabs(x - xsave) > fabs(ode_dxsave))
          if (ode_count < ode_max - 1)
            {
              ode_xp[++ode_count] = x;
              for (int i = 1; i <= nvar; ++i)
                ode_yp[i][ode_count] = y[i];
              xsave = x;
            }
      if ((x + h - x2) * (x + h - x1) > 0.0)
        h = x2 - x;
      (*stepper)(y, dydx, nvar, &x, h, eps, yscale, &hdid, &hnext, derivs);
      if (hdid == h)
        ++(*nok);
      else
        ++(*nbad);
      if ((x - x2) * (x2 - x1) >= 0.0)
        {
          for (int i = 1; i <= nvar; ++i)
            ystart[i] = y[i]; 
          if (ode_max)
            {
              ode_xp[++ode_count] = x;
              for (int i = 1; i <= nvar; ++i)
                ode_yp[i][ode_count] = y[i];
            }
          free_dvector(dydx, 1, nvar);
          free_dvector(y, 1, nvar);
          free_dvector(yscale, 1, nvar);
          return;
        }
      if (fabs(hnext) <= hmin)
        nrerror("Step size to small in ode_integrate.");
      h = hnext;
  }
  nrerror("Too many steps in routine ode_integrate.");
}



///
///  Modified midpoint step.  At xs, input the dependent variable vector y[1..nvar],
///  and its derivative dydx[1..nvar].  Also input is htot, the total step to be made,
///  and nstep, the number of interior steps to be used.  The output is returned as 
///  yout[1..nvar], which need not be distinct from y[1..nvar]; if it is distinct
///  however, then y and dydx will be returned undamaged.  Derivs is the user-supplied
///  routine for calculating the right-hand side derivative.
///
void
modified_midpoint(double *y, double *dydx, int nvar, double xs,
                  double htot, int nstep, double *yout,
                  void (*derivs)(double, double *, double *))
{
  double * ym = dvector(1, nvar);
  double * yn = dvector(1, nvar);

  double h = htot / nstep;
  for (int i = 1; i <= nvar; ++i)
    {
      ym[i] = y[i];
      yn[i] = y[i] + h*dydx[i];
    }
  double x = xs + h;
  (*derivs)(x, yn, yout);
  double h2 = 2.0 * h;
  for (int n = 2; n <= nstep; ++n)
    {
      for (int i = 1; i <= nvar; ++i)
        {
          double swap = ym[i] + h2 * yout[i];
          ym[i] = yn[i];
          yn[i] = swap;
        }
      x += h;
      (*derivs)(x, yn, yout);
    }
  for (int i = 1; i <= nvar; ++i)
    yout[i] = 0.5 * (ym[i] + yn[i] + h * yout[i]);
  free_dvector(yn, 1, nvar);
  free_dvector(ym, 1, nvar);
}


///
///  Stoermer's rule for integrating second order conservative systems of the form
///  y'' = f(x,y) for a system of n = nv/2 equations.  On input y[1..nv] contains
///  y in the first n elements and y' in the second n elements all evaluated at xs.
///  d2y[1..nv] contains the right hand side function f (also evaluated at xs) in
///  its first n elements (the second n elements are not referenced).  Also input
///  is htot, the total step to be taken and nstep, the number of substeps to be used.
///  The output is returned as yout[1..nv], with the same storage arrangement as y.
///  derivvs is the user-supplied routine that calculates f.
///
///  This routine can replace modified_midpoint above.
///
void
stoermer(double *y, double *d2y, int nv, double xs,
         double htot, int nstep, double *yout,
         void (*derivs)(double, double *, double *))
{
  double *ytemp = dvector(1, nv);

  double h = htot / nstep;
  double hh = 0.5 * h;
  int neqns = nv / 2;
  for (int i = 1; i <= neqns; ++i)
    {
      int n = neqns + i;
      ytemp[i] = y[i] + (ytemp[n] = h * (y[n] + hh * d2y[i]));
    }
  double x = xs + h;
  (*derivs)(x, ytemp, yout);
  double h2 = 2.0 * h;
  for (int nn = 2; nn <= nstep; ++nn)
    {
      for (int i = 1; i <= neqns; ++i)
        {
	  int n = neqns + i;
          ytemp[i] += (ytemp[n] += h2 * yout[i]);
	}
      x += h;
      (*derivs)(x, ytemp, yout);
    }
  for (int i = 1; i <= neqns; ++i)
    {
      int n = neqns + i;
      yout[n] = ytemp[n] / h + hh * yout[i];
      yout[i] = ytemp[i];
    }

  free_dvector(ytemp, 1, nv);
}




/*
 *    Global variables for the storage of extrapolation results for bulirsch_stoer.
 */
static double *bs_xx, **bs_yy;

static void bs_rat_extrap(int iest, double xest, double *yest,
                          double *yz, double *dy, int nv);

static void bs_poly_extrap(int iest, double xest, double *yest,
                           double *yz, double *dy, int nv);


#define KMAXX 9
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALEMAX 0.1

///
///  Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy
///  and adjust stepsize.  Input are the dependent variables y[1..nv] and the derivatives dydx[1..nv]
///  at the starting value of the independent variable x.
///
void
bulirsch_stoer(double *y, double *dydx, int nv, double *xx,
               double htry, double eps, double *yscale,
               double *hdid, double *hnext,
               void (*derivs)(double, double *, double *))
{
    int iq, k, kk, km, reduct, exitflag = 0;
    double eps1, errmax, fact, h, red, scale, work, workmin, xest, *err, *yerr, *ysav, *yseq;
    static int first = 1, kmax, kopt, nseq[IMAXX+1] = { 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20 };
    static double epsold = -1.0, xnew, a[IMAXX+1], alf[KMAXX+1][KMAXX+1];

    bs_xx = dvector(1, KMAXX);
    bs_yy = dmatrix(1, nv, 1, KMAXX);
    err = dvector(1, KMAXX);
    yerr = dvector(1, nv);
    ysav = dvector(1, nv);
    yseq = dvector(1, nv);

    if (eps != epsold)
      {
        *hnext = xnew = -1.0e29;
        eps1 = SAFE1*eps;
        a[1] = nseq[1] + 1;
        for (int k = 1; k <= KMAXX; ++k)
          a[k+1] = a[k] + nseq[k+1];
        for (int iq = 2; iq <= KMAXX; ++iq)
          for (int k = 1; k < iq; ++k)
            alf[k][iq] = pow(eps1, (a[k + 1] - a[iq + 1]) / ((a[iq+1] - a[1]) * (2 * k + 1)));
        epsold = eps;
        for (int kopt = 2; kopt < KMAXX; ++kopt)
          if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt])
            break;
        kmax = kopt;
      }
    h = htry;
    for (int i = 1; i <= nv; ++i)
      ysav[i] = y[i];
    if (*xx != xnew || h != *hnext)
      {
        first = 1;
        kopt = kmax;
      }
    reduct = 0;
    while (1)
      {
        for (int k = 1; k <= kmax; ++k)
          {
            xnew = *xx + h;
            if (xnew == *xx)
              nrerror("Step size underflow in bulirsch_stoer.");
            modified_midpoint(ysav, dydx, nv, *xx, h, nseq[k], yseq, derivs);
            xest = dsqr(h / nseq[k]);
            bs_poly_extrap(k, xest, yseq, y, yerr, nv);
            if (k != 1)
              {
                errmax = TINY;
                for (int i = 1; i <= nv; ++i)
                  errmax = dmax(errmax, fabs(yerr[i] / yscale[i]));
                errmax /= eps;
                km = kk - 1;
                err[km] = pow(errmax / SAFE1, 1.0 / (2 * k + 1));
              }
            if (k != 1 && (k >= kopt - 1 || first))
              {
                if (errmax < 1.0)
                  {
                    exitflag = 1;
                    break;
                  }
                if (k == kmax || k == kopt + 1)
                  {
                    red = SAFE2 / err[km];
                    break;
                  }
                else if (k == kopt && alf[kopt - 1][kopt] < err[km])
                  {
                    red = 1.0 / err[km];
                    break;
                  }
                else if (k == kopt && alf[km][kmax - 1] < err[km])
                  {
                    red = alf[km][kopt - 1] / err[km];
                    break;
                  }
              }
          }
        if (exitflag)
          break;
        red = dmin(red, REDMIN);
        red = dmax(red, REDMAX);
        h *= red;
        reduct = 1;
      }
    *xx = xnew;
    *hdid = h;
    first = 0;
    workmin = 1.0e35;
    for (int kk = 1; kk <= km; ++kk)
      {
        fact = dmax(err[kk], SCALEMAX);
        work = fact * a[kk + 1];
        if (work < workmin)
          {
            scale = fact;
            workmin = work;
            kopt = kk+1;
          }
      }
    *hnext = h / scale;
    if (kopt >= k && kopt != kmax && !reduct)
      {
        fact = dmax(scale/alf[kopt-1][kopt], SCALEMAX);
        if (fact * a[kopt + 1] <= workmin)
          {
            *hnext = h/fact;
            ++kopt;
          }
      }

    free_dvector(yseq, 1, nv);
    free_dvector(ysav, 1, nv);
    free_dvector(yerr, 1, nv);
    free_dvector(err, 1, KMAXX);
    free_dmatrix(bs_yy, 1, nv, 1, KMAXX);
    free_dvector(bs_xx, 1, KMAXX);
}




/*
 *    Routine used by bulirsch_stoer to perform rational function extrapolation.
 */
static void
bs_rat_extrap(int iest, double xest, double *yest, 
              double *yz, double *dy, int nv)
{
  double * fx = dvector(1, iest);

  bs_xx[iest] = xest;
  if (iest == 1)
    for (int j = 1; j <= nv; ++j)
      dy[j] = bs_yy[j][1] = yz[j] = yest[j];
  else
    {
      //  Evaluate next diagonal in the tableau.
      for (int k = 1; k < iest-1; ++k)
        fx[k + 1] = bs_xx[iest-k]/xest;
      for (int j = 1; j <= nv; ++j)
        {
          double v = bs_yy[j][1];
          double yy, c, ddy;
          bs_yy[j][1] = c = yy = yest[j];
          for (int k = 2; k < iest; ++k)
            {
              double b1 = fx[k] * v;
              double b = b1 - c;
              //  Watch division by zero.
              if (b)
        	{
                  b = (c - v) / b;
                  ddy = c * b;
                  c = b1 * b;
        	}
              else
        	ddy = v;
              if (k != iest)
        	v = bs_yy[j][k];
              bs_yy[j][k] = ddy;
              yy += ddy;
            }
          dy[j] = ddy;
          yz[j] = yy;
        }
    }
  free_dvector(fx, 1, iest);
}




/*
 *    Routine used by bulirsch_stoer to perform polynomial function extrapolation.
 */
static void
bs_poly_extrap(int iest, double xest, double *yest, 
               double *yz, double *dy, int nv)
{
  double * c = dvector(1, nv);

  bs_xx[iest] = xest;
  for (int j = 1; j <= nv; ++j)
    dy[j] = yz[j] = yest[j];
  if (iest == 1)
    for (int j = 1; j <= nv; ++j) bs_yy[j][1] = yest[j];
  else
    {
      for (int j = 1; j <= nv; ++j)
        c[j] = yest[j];
      for (int k = 1; k <= nv; ++k)
        {
          double delta = 1.0 / (bs_xx[iest-k] - xest);
          double f1 = delta * xest;
          double f2 = delta * bs_xx[iest - k];
          for (int j = 1; j <= nv; ++j)
            {
              double q = bs_yy[j][k];
              bs_yy[j][k] = dy[j];
              delta = c[j] - q;
              dy[j] = delta*f1;
              c[j] = delta*f2;
              yz[j] += dy[j];
            }
        }
      for (int j = 1; j <= nv; ++j)
        bs_yy[j][iest] = dy[j];
    }
  free_dvector(c, 1, nv);
}


