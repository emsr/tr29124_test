/*
g++ -std=c++14 -I. -o test_airy_fock test_airy_fock.cpp
./test_airy_fock > test_airy_fock.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>

#include <airy_fock.h>

int
main()
{
  using namespace std::complex_literals;
  constexpr double dpi = 3.141592653589793238462643383279502884195;
  using cmplx = std::complex<double>;

  std::cout.precision(std::numeric_limits<double>::digits10);
  const auto w = 6 + std::cout.precision();
  const int mm = 20;

  for (int ii = -mm; ii <= mm; ++ii)
    {
      if (ii > -mm)
        std::cout << '\n';
      for (int jj = -mm; jj <= mm; ++jj)
        {
          const auto t = cmplx(ii * 1.0, jj * 1.0);

          cmplx vt1, vt2;
          airy_fock(FockAiryType::Value, t, vt1, vt2);

          cmplx dt1, dt2;
          airy_fock(FockAiryType::Deriv, t, dt1, dt2);

          auto dW = -0.5i * (vt2 * dt1 - dt2 * vt1);

          std::cout << ' ' << std::setw(w) << t.real() << ' ' << std::setw(w) << t.imag()
                    << ' ' << std::setw(w) << vt1.real() << ' ' << std::setw(w) << vt1.imag() << ' ' << std::setw(w) << vt2.real() << ' ' << std::setw(w) << vt2.imag()
                    << ' ' << std::setw(w) << dt1.real() << ' ' << std::setw(w) << dt1.imag() << ' ' << std::setw(w) << dt2.real() << ' ' << std::setw(w) << dt2.imag()
                    << ' ' << std::setw(w) << dW.real() << ' ' << std::setw(w) << dW.imag()
                    << '\n';
        }
    }
}
