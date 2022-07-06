#include <iostream>
#include <iomanip>

#include "CohenSum.h"

template<typename Real>
  int
  test_CohenSum()
  {
    constexpr auto w = 1 + std::numeric_limits<Real>::max_digits10;
    std::cout.precision(w);

    emsr::CohenSum<Real> CS;

    constexpr auto pipid12 = Real{0.8224670334241132182362075833230125946103L};
    Real naive = 0, sign = 1;
    for (unsigned int k = 0; k < CS.max_num_terms(); ++k)
    {
      CS += Real{1} / ((k + 1) * (k + 1));
      naive += sign / ((k + 1) * (k + 1));
      sign = -sign;
    }
    std::cout << "  CS    = " << std::setw(w) << CS() << '\n';
    std::cout << "  naive = " << std::setw(w) << naive << '\n';
    std::cout << "  true  = " << std::setw(w) << pipid12 << '\n';
    std::cout << "  err   = " << std::setw(w) << CS() - pipid12 << '\n';
    return 0;
  }

int
main()
{
  int num_errors = 0;
  std::cout << "\n  float";
  std::cout << "\n  -----\n";
  num_errors += test_CohenSum<float>();

  std::cout << "\n  double";
  std::cout << "\n  ------\n";
  num_errors += test_CohenSum<double>();

  std::cout << "\n  long double";
  std::cout << "\n  -----------\n";
  num_errors += test_CohenSum<long double>();
}
