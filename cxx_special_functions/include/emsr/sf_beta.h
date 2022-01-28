#ifndef SF_BETA_H
#define SF_BETA_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_beta.tcc>

namespace emsr
{

  // Beta functions

  /**
   * Return the beta function, @f$ B(a,b) @f$, for @c float parameters
   * @c a, @c b.
   *
   * @see beta for more details.
   */
  inline float
  betaf(float a, float b)
  { return emsr::detail::beta<float>(a, b); }

  /**
   * Return the beta function, @f$B(a,b)@f$, for long double
   * parameters @c a, @c b.
   *
   * @see beta for more details.
   */
  inline long double
  betal(long double a, long double b)
  { return emsr::detail::beta<long double>(a, b); }

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
   * Return the regularized incomplete beta function of parameters
   * @c a, @c b, and argument @c x.
   *
   * See ibeta for details.
   */
  inline float
  ibetaf(float a, float b, float x)
  { return emsr::detail::beta_inc<float>(a, b, x); }

  /**
   * Return the regularized incomplete beta function of parameters @c a, @c b,
   * and argument @c x.
   *
   * See ibeta for details.
   */
  inline long double
  ibetal(long double a, long double b, long double x)
  { return emsr::detail::beta_inc<long double>(a, b, x); }

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

  inline float
  ibetacf(float a, float b, float x)
  { return 1.0F - emsr::ibetaf(a, b, x); }

  inline long double
  ibetacl(long double a, long double b, long double x)
  { return 1.0L - emsr::ibetal(a, b, x); }

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
