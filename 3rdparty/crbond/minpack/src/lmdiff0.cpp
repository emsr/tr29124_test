/*
 * lmdif0.cpp -- driver for lmdif
 */

#include <cmath>
#include <vector>
#include "cminpak.h"

int lmdif0(void fcn(), int m, int n, double x[], int msk[],
           double fvec[], double tol, int &info, int &nfev)
{
 // Check input parameters.
 if (n <= 0 || m < n || tol < 0.0)
   {
      info = 0;
      return 1;
   }
  // Allocate memory for working arrays.
  std::vector<int> ipvt(n);
  std::vector<double> diag(n);
  std::vector<double> qtf(n);
  std::vector<double> wa1(n);
  std::vector<double> wa2(n);
  std::vector<double> wa3(n);
  std::vector<double> wa4(m);


  // Create 2d matrix for Jacobian.
  std::vector<std::vector<double>> fjac(n, std::vector<double>(m));

  // Set convergence tolerances.
  auto ftol = tol;
  auto xtol = tol;
  auto gtol = 0.0;

  auto maxfev = 800;
  auto epsfcn = 0.0;
  auto mode = 1;
  auto factor = 10;
  auto nfev = 0;

  lmdif(fcn,m,n,x,msk,fvec,ftol,xtol,gtol,maxfev,epsfcn,diag,mode,
        factor,info,nfev,fjac,ipvt,qtf,wa1,wa2,wa3,wa4);

  if (info == 8)
    info = 4;

  return 0;
}
    



    
