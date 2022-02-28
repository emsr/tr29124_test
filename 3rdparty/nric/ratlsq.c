
#include <stdio.h>
#include <math.h>
#include "nric.h"

#define NPFAC 8
#define MAXIT 5
#define BIG 1.0e30

double
ratval(double x, double cof[], int mm, int kk)
{
    double sumn = cof[mm];
    for (int j = mm - 1; j >= 0; --j)
        sumn = sumn * x + cof[j];
    double sumd = 0.0;
    for (int j = mm + kk; j >= mm + 1; --j)
        sumd = (sumd + cof[j]) * x;
    return sumn / (1.0 + sumd);
}

void
ratlsq(double (*fn)(double), double a, double b, int mm, int kk, double cof[], double *dev)
{
    int ncof, npt;
    double devmax, e, hth, power, sum, *bb, *coff, *ee, *fs, **u, **v, *w, *wt, *xs;
    ncof = mm + kk + 1;
    npt = NPFAC * ncof;
    bb = dvector(1,npt);
    coff = dvector(0,ncof-1);
    ee = dvector(1,npt);
    fs = dvector(1,npt);
    u = dmatrix(1,npt,1,ncof);
    v = dmatrix(1,ncof,1,ncof);
    w = dvector(1,ncof);
    wt = dvector(1,npt);
    xs = dvector(1,npt);

    *dev = BIG;
    for (int i = 1; i <= npt; ++i)
    {
        if (i < npt/2)
        {
            hth = PIO2*(i-1)/(npt-1.0);
            xs[i] = a + (b - a) * dsqr(sin(hth));
        }
        else
        {
            hth = PIO2*(npt-i)/(npt-1.0);
            xs[i] = b - (b-a)*dsqr(sin(hth));
        }
        fs[i] = (fn)(xs[i]);
        wt[i] = 1.0;
        ee[i] = 1.0;
    }

    e = 0.0;
    for (int it = 1; it <= MAXIT; ++it)
    {
        for (int i = 1; i <= npt; ++i)
        {
            power = wt[i];
            bb[i] = power * (fs[i] + dsign(e, ee[i]));
            for (int j = 1; j <= mm + 1; ++j)
            {
                u[i][j] = power;
                power *= xs[i];
            }
            power = -bb[i];
            for (int j = mm + 2; j <= ncof; ++j)
            {
                power *= xs[i];
                u[i][j] = power;
            }
        }
/*
        for (int i = 1; i <= npt; ++i)
        {
            for (int j = 1; j <= ncof; ++j)
                printf(" %f", u[i][j]);
            printf("\n");
        }
*/
        sv_decomp(u, npt, ncof, w, v);
        sv_backsub(u, w, v, npt, ncof, bb, coff - 1);

        devmax = sum = 0.0;
        for (int j = 1; j <= npt; ++j)
        {
            ee[j] = ratval(xs[j], coff, mm, kk) - fs[j];
            wt[j] = fabs(ee[j]);
            sum += wt[j];
            if (wt[j] > devmax)
                devmax = wt[j];
        }
        e = sum / npt;
        if (devmax <= *dev)
        {
            *dev = devmax;
            for (int j = 0; j < ncof; ++j)
                cof[j] = coff[j];
        }
        printf("# ratlsq iterations = %2d; max error = %10.3e\n", it, devmax);
        printf("# numer: ");
        for (int i = 0; i < mm + 1; ++i)
            printf(" %14.6f", coff[i]);
        printf("\n");
        printf("# denom: %14.6f", 1.0);
        for (int i = mm + 1; i < mm + kk + 1; ++i)
            printf(" %14.6f", coff[i]);
        printf("\n");

        printf("\n\n");
        for (int i = 1; i <= 100; ++i)
        {
            double x = 0.01 * i * PI;
            double r = ratval(x, cof, mm, kk);
            double f = (*fn)(x);
            printf(" %12.6f %12.6f %12.6f %12.6f\n", x, r, f, r - f);
        }
    }

    free_dvector(bb, 1,npt);
    free_dvector(coff, 0,ncof-1);
    free_dvector(ee, 1,npt);
    free_dvector(fs, 1,npt);
    free_dmatrix(u, 1,npt,1,ncof);
    free_dmatrix(v, 1,ncof,1,ncof);
    free_dvector(w, 1,ncof);
    free_dvector(wt, 1,npt);
    free_dvector(xs, 1,npt);
}

double
fun(double x)
{
  return cos(x) / (1.0 - exp(x));
}

int
main()
{
    int mm = 4;
    int kk = 4;
    double cof[10];
    double dev;
    ratlsq(fun, 0.01, PI, mm, kk, cof, &dev);
    printf("# error: %14.6f\n", dev);
    printf("# numer: ");
    for (int i = 0; i < mm + 1; ++i)
        printf(" %14.6f", cof[i]);
    printf("\n");
    printf("# denom: %14.6f", 1.0);
    for (int i = mm + 1; i < mm + kk + 1; ++i)
        printf(" %14.6f", cof[i]);
    printf("\n");

    printf("\n\n");
    for (int i = 1; i <= 100; ++i)
    {
        double x = 0.01 * i * PI;
        double r = ratval(x, cof, mm, kk);
        double f = fun(x);
        printf(" %12.6f %12.6f %12.6f %12.6f\n", x, r, f, r - f);
    }
}
