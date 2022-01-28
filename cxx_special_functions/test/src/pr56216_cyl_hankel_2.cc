
#include <limits>
#include <complex>

#include <emsr/special_functions.h>

int
test01()
{
  const double inf = std::numeric_limits<double>::infinity();
  std::complex<double> h20 = emsr::cyl_hankel_2(0.0, 0.0);
  double h20r = std::real(h20);
  double h20i = std::imag(h20);
  std::complex<double> h21 = emsr::cyl_hankel_2(1.0, 0.0);
  double h21r = std::real(h21);
  double h21i = std::imag(h21);

  int num_errors = 0;
  if (h20r != 1.0)
    ++num_errors;
  if (h20i != inf)
    ++num_errors;
  if (h21r != 0.0)
    ++num_errors;
  if (h21i != inf)
    ++num_errors;
  return num_errors;
}

int
main()
{
  return test01();
}
