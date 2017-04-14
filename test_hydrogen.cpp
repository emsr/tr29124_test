
//
//  This tests both the Laurerre polynomials and the spherical legendre function.
//

#include <complex>
#include <cmath>
#include <tr1/cmath>

#include "hydrogen.tcc"

template<typename _Tp>
  _Tp
  __coulomb_norm(unsigned int __l, _Tp __eta)
  {
    _Tp _Ck = std::sqrt(_S_2pi * __eta / (std::exp(_S_2pi * __eta) - _Tp{1}));
    if (__l == 0)
      return _Ck;
    else
      {
	for (int __k = 0; __k < __l; ++__k)
	  _Ck *= std::hypot(_Tp(__k + 1), __eta) / _Tp(__k) / _Tp(__k + 1);
	return _Ck;
      }
  }

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

