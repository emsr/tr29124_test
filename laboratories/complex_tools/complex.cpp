/**
 *
 */

#include <iostream>
#include <complex>
#include <tr1/cmath>

int
main()
{
  double nu = 0.0;
  std::complex<double> z(1.0, 1.0);

  //std::complex<double> w = std::tr1::cyl_bessel_j( 0.0, z );
  //std::cout << " z = " << z << "  w = " << w << '\n';

  std::complex<double> Jnu, Nnu, Jpnu, Npnu;
  emsr::detail::bessel_jn(std::complex<double>(nu), z, Jnu, Nnu, Jpnu, Npnu);
  std::cout << " z = " << z << '\n';
  std::cout << " Jnu = " << Jnu << '\n';
  std::cout << " Nnu = " << Nnu << '\n';
  std::cout << " Jpnu = " << Jpnu << '\n';
  std::cout << " Npnu = " << Npnu << '\n';

  return 0;
}
