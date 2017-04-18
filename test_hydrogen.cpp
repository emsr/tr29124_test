/*
$HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_hydrogen test_hydrogen.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_hydrogen > test_hydrogen.txt

$HOME/bin/bin/g++ -std=c++17 -std=gnu++17 -g -Wall -Wextra -I. -o test_hydrogen test_hydrogen.cpp -lquadmath -Lwrappers/debug -lwrap_gsl
PATH=wrappers/debug:$PATH ./test_hydrogen > test_hydrogen.txt
*/

//
//  This tests both the Laurerre polynomials and the spherical legendre function.
//

#include <complex>
#include <cmath>

template <typename _Tp>
  void
  test01()
  {
    const _Tp c = 2.99792458e8L;  //  Speed of light in m/s.
    const _Tp hbar = 1.0546e-27L;  //  Planck's constant in J s.
    const _Tp hbarc = hbar * c;
    const _Tp epsilon0 = 8.8542e-12L;  //  Permittivity of free space C^2 / J m.
    const _Tp mu0 = _Tp(4.0e-7L) * _Tp(__PI);  //  Permeability of free space in N / A^2.
    const _Tp e = 1.6e-19L;  //  Charge of the electron in C.
    const _Tp me = 9.1095e-31L;  // Mass of the electron in kg.
    const _Tp alpha = e * e / (_Tp(4) * _Tp(__PI) * epsilon0);
    const _Tp a0 = hbarc * hbarc / (me * c * c * alpha);  //  Bohr radius.

    return;
  }

template <typename _Tp>
  void
  test_norm()
  {
    for (unsigned int l = 0; l <= 100; ++l)
      for (int k = 0; k <= 20; ++k)
	{
	  auto eta = k * _Tp{1};
	  auto C_cxx = std::__detail::__coulomb_norm(l, eta);
	  auto C_gsl = gsl::coulomb_norm(l, eta);
	  std::cout << setw(4) << l
		    << setw(w) << eta
		    << setw(w) << C_cxx
		    << setw(w) << C_gsl
		    << '\n';
	}
  }

int
main()
{
}


