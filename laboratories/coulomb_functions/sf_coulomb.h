#ifndef SF_COULOMB_H
#define SF_COULOMB_H 1

namespace emsr
{
  template<typename Tp>
    struct Coulomb
    {
      Tp F_value;
      Tp F_deriv;
      Tp F_exp;
      Tp G_value;
      Tp G_deriv;
      Tp G_exp;
    };
}

#include "sf_coulomb.tcc"

namespace emsr
{
  template<typename Tp>
    Tp
    coulomb_norm(Tp lambda, Tp eta)
    {
      return emsr::detail::coulomb_norm(lambda, eta);
    }

  template<typename Tp>
    Tp
    coulomb_f(Tp lambda, Tp eta, Tp rho)
    {
      return emsr::detail::coulomb_wave_FG(lambda, eta, rho, 0).F_value;
    }

  template<typename Tp>
    Tp
    coulomb_g(Tp lambda, Tp eta, Tp rho)
    {
      return emsr::detail::coulomb_wave_FG(lambda, eta, rho, 0).G_value;
    }

} // namespace emsr

#endif // SF_COULOMB_H
