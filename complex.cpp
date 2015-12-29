
// $HOME/bin_specfun/bin/g++ -o complex complex.cpp

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./complex

#include <iostream>
#include <complex>
#include <tr1/cmath>

int
main()
{
  double nu = 0.0;
  std::complex<double> z(1.0, 1.0);

  //std::complex<double> w = std::tr1::cyl_bessel_j( 0.0, z );
  //std::cout << " z = " << z << "  w = " << w << std::endl;

  std::complex<double> Jnu, Nnu, Jpnu, Npnu;
  std::tr1::__detail::__bessel_jn(std::complex<double>(nu), z, Jnu, Nnu, Jpnu, Npnu);
  std::cout << " z = " << z << std::endl;
  std::cout << " Jnu = " << Jnu << std::endl;
  std::cout << " Nnu = " << Nnu << std::endl;
  std::cout << " Jpnu = " << Jpnu << std::endl;
  std::cout << " Npnu = " << Npnu << std::endl;

  return 0;
}
