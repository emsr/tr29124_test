
// Crash of Bessel functions at x==0!

#include <limits>
#include <complex>

#include <emsr/special_functions.h>

int
test01()
{
  const double inf = std::numeric_limits<double>::infinity();
  std::complex<double> h10 = emsr::cyl_hankel_1(0.0, 0.0);
  double h10r = std::real(h10);
  double h10i = std::imag(h10);
  std::complex<double> h11 = emsr::cyl_hankel_1(1.0, 0.0);
  double h11r = std::real(h11);
  double h11i = std::imag(h11);

  int num_errors = 0;
  if (h10r != 1.0)
    ++num_errors;
  if (h10i != -inf)
    ++num_errors;
  if (h11r != 0.0)
    ++num_errors;
  if (h11i != -inf)
    ++num_errors;
  return num_errors;
}

int
main()
{
  return test01();
}
