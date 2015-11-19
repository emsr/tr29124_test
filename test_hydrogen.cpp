
//
//  This tests both the Laurerre polynomials and the spherical legendre function.
//

#include <complex>
#include <cmath>
#include <tr1/cmath>

#include "hydrogen.tcc"

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

