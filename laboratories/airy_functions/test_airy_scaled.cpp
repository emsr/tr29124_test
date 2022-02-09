#include <iostream>

#include <emsr/sf_mod_bessel.h>
#include <wrap_gsl.h>

int
main()
{
  for (int i = 0; i <= 100; ++i)
  {
    auto x = 0.1 * i;
    auto xi = (2.0/3.0) * std::pow(x, 3.0/2.0);
    auto expxi = std::exp(xi);
    auto y = emsr::detail::airy(x, false);
    auto ys = emsr::detail::airy(x, true);
    auto Ai = gsl::airy_ai_scaled(x);
    auto Bi = gsl::airy_bi_scaled(x);
    std::cout << ' ' << x
              << ' ' << Ai
              << ' ' << y.Ai_value * expxi
              << ' ' << ys.Ai_value
              << ' ' << Bi
              << ' ' << y.Bi_value / expxi
              << ' ' << ys.Bi_value
              << '\n';
  }
}
