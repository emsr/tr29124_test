#include <iostream>
#include <iomanip>
#include <e_float/e_float.h>
#include <functions/functions.h>

int main(void)
{
  // Calculate some real numbered function values...
  static const e_float x = ef::pi() / 7;
  static const e_float s = ef::sin(x);
  static const e_float g = ef::gamma(x);
  static const e_float b = ef::cyl_bessel_j(ef::third(), x);

  // ...and print them to full precision.
  std::cout.precision(std::numeric_limits<e_float>::digits10);

  std::cout << s << std::endl;
  std::cout << g << std::endl;
  std::cout << b << std::endl;
}
