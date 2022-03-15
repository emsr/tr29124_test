
// assoc_legendre returns negated value when m is odd

#include <emsr/special_functions.h>

int
test01()
{
  long double P11_ok = 0.866025403784438646787;

  int num_errors = 0;
  auto P11 = emsr::assoc_legendre(1, 1, 0.5);
  auto diff_ok = P11 - P11_ok;
  if (std::abs(diff_ok) > 1.0e-12) ++num_errors;
  return num_errors;
}

int
main()
{
  return test01();
}
