#ifndef SF_STRUVE_TCC
#define SF_STRUVE_TCC 1

#include <emsr/math_constants.h>
#include <emsr/summation.h>
#include <emsr/special_functions.h>

namespace emsr
{
namespace detail
{

  /**
   * An enum to dispatch Struve function summation.
   */
  enum _StruveType
  : int
  {
    _StruveH,
    _StruveK,
    _StruveL,
    _StruveM
  };

  /**
   * Return either the Struve function of the first kind
   * @f$ \boldmath{H}_\nu(x) @f$ or the modified Struve function
   * of the first kind @f$ \boldmath{L}_\nu(x) @f$
   * depending on whether @c sign is -1 or +1 respectively.
   */
  template<_StruveType _Type, typename _Tp>
    _Tp
    struve_series(_Tp nu, _Tp x)
    {
      using _Val = _Tp;

      using _BasicSum = emsr::BasicSum<_Val>;
      using _WenigerBasSum = emsr::WenigerDeltaSum<_BasicSum>;
      using _WijnSum = emsr::VanWijngaardenSum<_Val>;
      using _WenigerWijnSum = emsr::WenigerDeltaSum<_WijnSum>;
      using _WenigerSum = std::conditional_t<_Type == _StruveH,
					     _WenigerWijnSum, _WenigerBasSum>;
      int sign = (_Type == _StruveH ? -1 : _Type == _StruveL ? +1 : 0);
      assert(sign != 0);

      constexpr int s_max_iter = 1000;
      const auto s_eps = emsr::epsilon(std::real(x));
      const auto s_sqrt_pi = emsr::sqrtpi_v<_Tp>;

      auto x2 = x / _Val{2};
      auto xx4 = _Tp(sign) * x2 * x2;
      auto term = _Val{1};
      auto struve = _WenigerSum(_Val{1});
      struve += term;
      for (int k = 1; k < s_max_iter; ++k)
	{
      	  term *= xx4 / _Val(k + 0.5L) / (nu + _Val(k + 0.5L));
	  struve += term;
	  if (std::abs(term) < s_eps * std::abs(struve))
	    break;
	}
      auto factor = _Val{2} * std::pow(x2, nu + _Val{1})
		    / emsr::tgamma(nu + _Val{1.5L}) / s_sqrt_pi;

      return factor * struve();
    }

  /**
   * Return either the Struve function of the second kind
   * @f$ \boldmath{K}_\nu(x) @f$ or the modified Struve function
   * of the second kind @f$ \boldmath{M}_\nu(x) @f$
   * depending on whether @c sign is +1 or -1 respectively.
   */
  template<_StruveType _Type, typename _Tp>
    _Tp
    struve_asymp(_Tp nu, _Tp x)
    {
      using _Val = _Tp;

      using _BasicSum = emsr::BasicSum<_Val>;
      using _WenigerBasSum = emsr::WenigerDeltaSum<_BasicSum>;
      using _WijnSum = emsr::VanWijngaardenSum<_Val>;
      using _WenigerWijnSum = emsr::WenigerDeltaSum<_WijnSum>;
      using _WenigerSum = std::conditional_t<_Type == _StruveM,
					     _WenigerWijnSum, _WenigerBasSum>;

      int sign = (_Type == _StruveK ? +1 : _Type == _StruveM ? -1 : 0);
      assert(sign != 0);

      constexpr int s_max_iter = 1000;
      const auto s_eps = emsr::epsilon(std::real(x));
      const auto s_sqrt_pi = emsr::sqrtpi_v<_Tp>;

      auto x2 = x / _Val{2};
      auto xx4 = _Val(sign) * x2 * x2;
      auto term = _Val{1};
      auto struve = _WenigerSum(_Val{1});
      struve += term;
      for (int k = 1; k < s_max_iter; ++k)
	{
	  const auto term_prev = term;
      	  term *= _Val(k - 0.5L) / (_Val(-k - 0.5L) + nu) / xx4;
	  if (std::abs(term) > std::abs(term_prev))
	    break;
	  struve += term;
	  if (std::abs(term) < s_eps * std::abs(struve))
	    break;
	}
      auto fact = _Val(sign) * std::pow(x2, nu - _Val{1})
		  / emsr::tgamma(nu + _Val{0.5L}) / s_sqrt_pi;

      return fact * struve();
    }

  /**
   * Return the Struve function of the first kind
   * @f$ \boldmath{H}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    struve_h(_Tp nu, _Tp x)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto s_nan = emsr::quiet_NaN(std::real(x));
      const auto s_max = emsr::digits10(std::real(x));

      if (std::real(x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	throw std::domain_error("struve_h: bad argument");
      else if (std::isnan(nu) || std::isnan(x))
	return s_nan;
      else if (std::abs(x) < s_max)
	return struve_series<_StruveH>(nu, x);
      else
	{
	  auto _Nnu = cyl_neumann_n(nu, x);
	  auto _Knu = struve_asymp<_StruveK>(nu, x);
	  return _Knu + _Nnu;
	}
    }

  /**
   * Return the Struve function of the second kind
   * @f$ \boldmath{K}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    struve_k(_Tp nu, _Tp x)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto s_nan = emsr::quiet_NaN(std::real(x));
      const auto s_max = emsr::digits10(std::real(x));

      if (std::real(x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	throw std::domain_error("struve_k: bad argument");
      else if (std::isnan(nu) || std::isnan(x))
	return s_nan;
      else if (std::abs(x) >= s_max)
	return struve_asymp<_StruveK>(nu, x);
      else
	{
	  auto _Nnu = cyl_neumann_n(nu, x);
	  auto _Hnu = struve_series<_StruveH>(nu, x);
	  return _Hnu - _Nnu;
	}
    }

  /**
   * Return the modified Struve function of the first kind
   * @f$ \boldmath{L}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    struve_l(_Tp nu, _Tp x)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto s_nan = emsr::quiet_NaN(std::real(x));
      const auto s_max = emsr::digits10(std::real(x));

      if (std::real(x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	throw std::domain_error("struve_l: bad argument");
      else if (std::isnan(nu) || std::isnan(x))
	return s_nan;
      else if (std::abs(x) < s_max)
	return struve_series<_StruveL>(nu, x);
      else
	{
	  auto _Inu = cyl_bessel_i(nu, x);
	  auto _Mnu = struve_asymp<_StruveM>(nu, x);
	  return _Mnu + _Inu;
	}
    }

  /**
   * Return the modified Struve function of the second kind
   * @f$ \boldmath{M}_\nu(x) @f$.
   */
  template<typename _Tp>
    _Tp
    struve_m(_Tp nu, _Tp x)
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto s_nan = emsr::quiet_NaN(std::real(x));
      const auto s_max = emsr::digits10(std::real(x));

      if (std::real(x) < _Real{0}) /// @todo Find out about Struve for x < 0.
	throw std::domain_error("struve_k: bad argument");
      else if (std::isnan(nu) || std::isnan(x))
	return s_nan;
      else if (std::abs(x) >= s_max)
	return struve_asymp<_StruveM>(nu, x);
      else
	{
	  auto _Inu = cyl_bessel_i(nu, x);
	  auto _Lnu = struve_series<_StruveL>(nu, x);
	  return _Lnu - _Inu;
	}
    }

} // namespace detail
} // namespace emsr

#endif // SF_STRUVE_TCC
