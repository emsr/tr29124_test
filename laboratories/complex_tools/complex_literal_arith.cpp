/**
 *
 */

#include <complex>

std::complex<double>
func(std::complex<double> z)
{
  return 2.0 * z;
}

void
test01()
{
  using namespace std::complex_literals;

  std::complex<double> z, nu;

  nu = 5.0;
  z = 1.0 - 3.0i;

  constexpr std::complex<double> H1{1.0 - 2.0i}, H1p{5.0 - 6.0i},
				 H2{3.0 - 4.0i}, H2p{7.0 - 8.0i};

  auto w = z - 2.0i;
  auto mu = 4.0i / nu;

  auto J = (H1 + H2) / 2.0;
  auto Jp = (H1p + H2p) / 2.0;
  auto Y = (H1 - H2) / 2.0i;
  auto Yp = (H1p - H2p) / 2.0i;

  auto Z = 1.5 * H1 + 3.5i * H2;

  auto f1 = func(z);

  auto f2 = func(3.5i);
}

int
main()
{
  test01();
}
