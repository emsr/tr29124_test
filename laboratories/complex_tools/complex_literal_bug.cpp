/**
 *
 */

#include <complex>

template<typename Real>
  std::complex<Real>
  func(std::complex<Real> z)
  {
    return Real{2} * z;
  }

void
test_float()
{
  using namespace std::literals::complex_literals;
  auto w = func(1.0f + 3.0if);
  w /= 2.0if;
  auto k = w / 4.0if;
  auto q = w * 1.5if;
}

void
test_double()
{
  using namespace std::literals::complex_literals;
  auto w = func(1.0 + 3.0i);
  w /= 2.0i;
  auto k = w / 4.0i;
  auto q = w * 1.5i;
}

void
test_long_double()
{
  using namespace std::literals::complex_literals;
  auto w = func(1.0l + 3.0il);
  w /= 2.0il;
  auto k = w / 4.0il;
  auto q = w * 1.5il;
}

int
main()
{
  test_float();
  test_double();
  test_long_double();
}
