/**
 *
 */

//
//  This tests both the Laurerre polynomials and the spherical legendre function.
//

#include <complex>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <emsr/math_constants.h>
#include <emsr/fp_type_util.h>
#include <emsr/complex_util.h>
#include <emsr/numeric_limits.h>
#include <emsr/specfun_state.h>
#include <emsr/sf_coulomb.h>
#include <wrap_gsl.h>

template <typename Tp>
  void
  test01()
  {
    const auto s_pi = emsr::pi_v<Tp>;
    const Tp c = Tp{2.99792458e8L};  //  Speed of light in m/s.
    const Tp hbar = Tp{1.0546e-27L};  //  Planck's constant in J s.
    const Tp hbarc = hbar * c;
    const Tp epsilon0 = Tp{8.8542e-12L};  //  Permittivity of free space C^2 / J m.
    const Tp mu0 = Tp{4.0e-7L} * s_pi;  //  Permeability of free space in N / A^2.
    const Tp e = Tp{1.6e-19L};  //  Charge of the electron in C.
    const Tp me = Tp{9.1095e-31L};  // Mass of the electron in kg.
    const Tp alpha = e * e / (Tp{4} * s_pi * epsilon0);
    const Tp a0 = hbarc * hbarc / (me * c * c * alpha);  //  Bohr radius.

    return;
  }

template <typename Tp>
  void
  test_norm()
  {
    std::cout.precision(8);
    std::cout.flags(std::ios::showpoint);
    const auto w = 8 + std::cout.precision();

    for (unsigned int l = 0; l <= 100; ++l)
      for (int k = 0; k <= 20; ++k)
	{
	  auto eta = k * Tp{1};
	  auto C_cxx = emsr::coulomb_norm(l, eta);
	  auto C_gsl = gsl::coulomb_norm(l, eta);
	  std::cout << std::setw(4) << l
		    << std::setw(w) << eta
		    << std::setw(w) << C_cxx
		    << std::setw(w) << C_gsl
		    << '\n';
	}
  }

int
main()
{
}


