#include "scorer.h"
#include <iostream>
#include <iomanip>

void
test_scorer()
{
  std::cout.precision(16);
  auto w = 6 + std::cout.precision();
  int ifac = 1;

  std::cout << "\n Scorer Gi(z)\n";
  for (int ix = 1; ix < 10; ++ix)
    for (int iy = 1; iy < 10; ++iy)
      {
	double x = 0.5 * ix;
	double y = 0.5 * iy;
	std::cout << " (" << x << ',' << y << ')';
	int ierr;
	double regi, imgi, regip, imgip;	
	scorer_gi(ifac, x, &y, &regi, &imgi, &regip, &imgip, &ierr);
	std::cout << " (" << std::setw(w) << regi << ',' << std::setw(w) << imgi << ')';
	std::cout << " (" << std::setw(w) << regip << ',' << std::setw(w) << imgip << ')';
	std::cout << '\n';
      }

  std::cout << "\n Scorer Hi(z)\n";
  for (int ix = 1; ix < 10; ++ix)
    for (int iy = 1; iy < 10; ++iy)
      {
	double x = 0.5 * ix;
	double y = 0.5 * iy;
	std::cout << " (" << x << ',' << y << ')';
	int ierr;
	double rehi, imhi, rehip, imhip;	
	scorer_hi(ifac, x, &y, &rehi, &imhi, &rehip, &imhip, &ierr);
	std::cout << " (" << std::setw(w) << rehi << ',' << std::setw(w) << imhi << ')';
	std::cout << " (" << std::setw(w) << rehip << ',' << std::setw(w) << imhip << ')';
	std::cout << '\n';
      }
}

int
main()
{
  test_scorer();
}
