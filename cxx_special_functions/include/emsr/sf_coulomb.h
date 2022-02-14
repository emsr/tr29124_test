#ifndef SF_COULOMB_H
#define SF_COULOMB_H 1

#include <complex>

#include <emsr/fp_type_util.h>

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

#include <emsr/detail/sf_coulomb.tcc>

namespace emsr
{
  template<typename Tpl, typename Tph>
    emsr::fp_promote_t<Tpl, Tph>
    coulomb_norm(Tpl lambda, Tph eta)
    {
      using type = emsr::fp_promote_t<Tpl, Tph>;
      return emsr::detail::coulomb_norm<type>(lambda, eta);
    }

  /**
   * Return the regular Coulomb function
   * @f[
   *    F_\lambda(\eta, \rho)
   * @f]
   */
  template<typename Tpl, typename Tph, typename Tp>
    emsr::fp_promote_t<Tpl, Tph, Tp>
    coulomb_f(Tpl lambda, Tph eta, Tp rho)
    {
      using type = emsr::fp_promote_t<Tpl, Tph, Tp>;
      return emsr::detail::coulomb_wave_FG<type>(lambda, eta, rho, 0).F_value;
    }

  /**
   * Return the irregular Coulomb function
   * @f[
   *    G_\lambda(\eta, \rho)
   * @f]
   */
  template<typename Tpl, typename Tph, typename Tp>
    emsr::fp_promote_t<Tpl, Tph, Tp>
    coulomb_g(Tpl lambda, Tph eta, Tp rho)
    {
      using type = emsr::fp_promote_t<Tpl, Tph, Tp>;
      return emsr::detail::coulomb_wave_FG<type>(lambda, eta, rho, 0).G_value;
    }

  /**
   * Return the Coulomb function
   * @f[
   *    H_\lambda^{(+)}(\eta, \rho) = G_\lambda(\eta, \rho) + i F_\lambda(\eta, \rho)
   * @f]
   */
  template<typename Tpl, typename Tph, typename Tp>
    std::complex<emsr::fp_promote_t<Tpl, Tph, Tp>>
    coulomb_hp(Tpl lambda, Tph eta, Tp rho)
    {
      using type = emsr::fp_promote_t<Tpl, Tph, Tp>;
      const auto s_i = std::complex<type>{0, 1};
      const auto FG = emsr::detail::coulomb_wave_FG<type>(lambda, eta, rho, 0);
      return FG.G_value + s_i * FG.F_value;
    }

  /**
   * Return the Coulomb function
   * @f[
   *    H_\lambda^{(-)}(\eta, \rho) = G_\lambda(\eta, \rho) - i F_\lambda(\eta, \rho)
   * @f]
   */
  template<typename Tpl, typename Tph, typename Tp>
    std::complex<emsr::fp_promote_t<Tpl, Tph, Tp>>
    coulomb_hm(Tpl lambda, Tph eta, Tp rho)
    {
      using type = emsr::fp_promote_t<Tpl, Tph, Tp>;
      const auto s_i = std::complex<type>{0, 1};
      const auto FG = emsr::detail::coulomb_wave_FG<type>(lambda, eta, rho, 0);
      return FG.G_value - s_i * FG.F_value;
    }

} // namespace emsr

#endif // SF_COULOMB_H
