#include <iostream>
#include <iomanip>
#include <functions/complex/e_float_complex.h>
#include <functions/functions.h>

int main(void)
{
  const ef_complex z1(ef::quarter(), ef::third());
  const ef_complex z2 = efz::sqrt(z1) + ef::tenth();
  const ef_complex z3 = efz::riemann_zeta(z2);

  std::cout << std::setprecision(std::numeric_limits<e_float>::digits10)
            << z3
            << std::endl;

  std::cout << std::setprecision(std::numeric_limits<e_float>::digits10)
            << efz::abs(z3)
            << std::endl;
}
