#ifndef SF_BETA_H
#define SF_BETA_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_beta.tcc>

namespace emsr
{

  // Beta functions

  /**
   * Return the beta function, @f$B(a,b)@f$, for real parameters
   * @c a, @c b.
   *
   * The beta function is defined by
   * @f[
   *   B(a,b) = \int_0^1 t^{a - 1} (1 - t)^{b - 1} dt
   *          = \frac{\Gamma(a)\Gamma(b)}{\Gamma(a+b)}
   * @f]
   * where @f$ a > 0 @f$ and @f$ b > 0 @f$
   *
   * @tparam _Tpa The floating-point type of the parameter @c a.
   * @tparam _Tpb The floating-point type of the parameter @c b.
   * @param a The first argument of the beta function, <tt> a > 0 </tt>.
   * @param b The second argument of the beta function, <tt> b > 0 </tt>.
   * @throw std::domain_error if <tt> a < 0 </tt> or <tt> b < 0 </tt>.
   */
  template<typename _Tpa, typename _Tpb>
    inline emsr::fp_promote_t<_Tpa, _Tpb>
    beta(_Tpa a, _Tpb b)
    {
      using type = emsr::fp_promote_t<_Tpa, _Tpb>;
      return emsr::detail::beta<type>(a, b);
    }

  // Incomplete beta functions

  /**
   * Return the regularized incomplete beta function of parameters @c a, @c b,
   * and argument @c x.
   *
   * The regularized incomplete beta function is defined by
   * @f[
   *    I_x(a, b) = \frac{B_x(a,b)}{B(a,b)}
   * @f]
   * where
   * @f[
   *   B_x(a,b) = \int_0^x t^{a - 1} (1 - t)^{b - 1} dt
   * @f]
   * is the non-regularized incomplete beta function and @f$ B(a,b) @f$
   * is the usual beta function.
   *
   * @param a The first parameter
   * @param b The second parameter
   * @param x The argument
   */
  template<typename Ta, typename Tb, typename Tp>
    inline emsr::fp_promote_t<Ta, Tb, Tp>
    ibeta(Ta a, Tb b, Tp x)
    {
      using type = emsr::fp_promote_t<Ta, Tb, Tp>;
      return emsr::detail::beta_inc<type>(a, b, x);
    }

  // Complementary incomplete beta functions

  /**
   * Return the regularized complementary incomplete beta function
   * of parameters @c a, @c b, and argument @c x.
   *
   * The regularized complementary incomplete beta function is defined by
   * @f[
   *    I_x(a, b) = I_x(a, b)
   * @f]
   *
   * @param a The parameter
   * @param b The parameter
   * @param x The argument
   */
  template<typename Ta, typename Tb, typename Tp>
    inline emsr::fp_promote_t<Ta, Tb, Tp>
    ibetac(Ta a, Tb b, Tp x)
    {
      using type = emsr::fp_promote_t<Ta, Tb, Tp>;
      return type(1) - emsr::ibeta<type>(a, b, x);
    }

} // namespace emsr

#endif // SF_BETA_H
