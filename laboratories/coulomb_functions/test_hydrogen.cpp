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
#include <ext/math_const.h>
#include <bits/numeric_limits.h>
#include <bits/specfun_state.h>
#include <bits/sf_trig.tcc>
#include <bits/sf_bernoulli.tcc>
#include <bits/sf_gamma.tcc>
#include <bits/sf_coulomb.tcc>
#include <wrap_gsl.h>

template <typename _Tp>
  void
  test01()
  {
    const auto _S_pi = __gnu_cxx::math::__pi_v<_Tp>;
    const _Tp c = _Tp{2.99792458e8L};  //  Speed of light in m/s.
    const _Tp hbar = _Tp{1.0546e-27L};  //  Planck's constant in J s.
    const _Tp hbarc = hbar * c;
    const _Tp epsilon0 = _Tp{8.8542e-12L};  //  Permittivity of free space C^2 / J m.
    const _Tp mu0 = _Tp{4.0e-7L} * _S_pi;  //  Permeability of free space in N / A^2.
    const _Tp e = _Tp{1.6e-19L};  //  Charge of the electron in C.
    const _Tp me = _Tp{9.1095e-31L};  // Mass of the electron in kg.
    const _Tp alpha = e * e / (_Tp{4} * _S_pi * epsilon0);
    const _Tp a0 = hbarc * hbarc / (me * c * c * alpha);  //  Bohr radius.

    return;
  }

template <typename _Tp>
  void
  test_norm()
  {
    std::cout.precision(8);
    std::cout.flags(std::ios::showpoint);
    const auto w = 8 + std::cout.precision();

    for (unsigned int l = 0; l <= 100; ++l)
      for (int k = 0; k <= 20; ++k)
	{
	  auto eta = k * _Tp{1};
	  auto C_cxx = std::__detail::__coulomb_norm(l, eta);
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


