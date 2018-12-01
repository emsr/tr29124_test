/*
g++ -std=c++2a -o run_coulfg coulfg.cpp run_coulfg.cpp
*/

struct steed_info
{
  double paccq = 0.0;
  int nfp = 0;
  int npq = 0;
  int iexp = 0;
  int m1 = 0;
};

void
coulfg(double xx,double eta1,double xlmin,double xlmax,
       double (&fc)[100], double (&gc)[100], double (&fcp)[100], double (&gcp)[100],
       int mode1, int kfn, int& ifail, steed_info& steed);

#include <iostream>
#include <iomanip>

int
main()
{
  std::cout << std::showpoint;

  int ifail;
  double fc[100], gc[100], fcp[100], gcp[100];

  steed_info steed;

  int mode = 1;
  int kfn = 0;
  double lam = 0.0;
  for (int ih = -2; ih <= 10; ih += 2)
    {
      auto eta = ih * 1.0;
      std::cout << "\n\neta = " << eta << '\n';
      for (int ir = 1; ir <= 200; ++ir)
	{
	  auto rho = ir * 0.1;
	  fc[0] = 0.0;
	  gc[0] = 0.0;
	  fcp[0] = 0.0;
	  gcp[0] = 0.0;
	  coulfg(rho, eta, lam, lam + 1.0, fc, gc, fcp, gcp, mode, kfn, ifail, steed);
	  std::cout << std::setw(16) << rho
                    << std::setw(16) << fc[1]
                    << std::setw(16) << gc[1]
                    << std::setw(16) << fcp[1]
                    << std::setw(16) << gcp[1]
                    << std::setw(4) << ifail
                    << std::setw(4) << steed.m1
                    << '\n';
	}
    }
}
