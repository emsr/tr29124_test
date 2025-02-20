/**
 *
 */

#include <cmath>
#include <complex>

#include <emsr/sf_legendre.h>

int
main()
{
  emsr::legendre(1, std::complex<double>{});
  emsr::legendre_q(1, std::complex<double>{});
  emsr::assoc_legendre(1, 0, std::complex<double>{});
  emsr::assoc_legendre_q(1, 0, std::complex<double>{});
}
