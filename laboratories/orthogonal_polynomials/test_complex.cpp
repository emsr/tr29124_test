/**
 *
 */

#include <cmath>
#include <complex>

int
main()
{
  std::legendre(1, std::complex<double>{});
  __gnu_cxx::legendre_q(1, std::complex<double>{});
  std::assoc_legendre(1, 0, std::complex<double>{});
  __gnu_cxx::assoc_legendre_q(1, 0, std::complex<double>{});
}
