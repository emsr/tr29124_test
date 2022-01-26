/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

//namespace emsr
//{
//namespace detail
//{

  /**
   * @brief  Return the Tricomi confluent hypergeometric function
   * @f[
   *   U(a,c,x) = \frac{\Gamma(1-c)}{\Gamma(a-c+1)} {}_1F_1(a;c;x)
   *       + \frac{\Gamma(c-1)}{\Gamma(a)} x^{1-c} {}_1F_1(a-c+1;2-c;x)
   * @f]
   * @param  a  The @a numerator parameter.
   * @param  c  The @a denominator parameter.
   * @param  x  The argument of the confluent hypergeometric function.
   * @return  The Tricomi confluent hypergeometric function.
   */
  template<typename _Tp>
    _Tp
    tricomi_u_naive(_Tp a, _Tp c, _Tp x)
    {
      auto U1 = _Tp{};
      auto b = a - c + _Tp{1};
      auto ib = emsr::fp_is_integer(b);
      if (!ib || ib() > 0)
	U1 = emsr::detail::gamma(_Tp{1} - c)
	     * emsr::detail::conf_hyperg(a, c, x)
	     / emsr::detail::gamma(b);

      auto U2 = _Tp{};
      auto ia = emsr::fp_is_integer(a);
      if (!ia || ia() > 0)
	U2 = emsr::detail::gamma(c - _Tp{1})
	     * std::pow(x, _Tp{1} - c)
	     * emsr::detail::conf_hyperg(b, _Tp{2} - c, x)
	     / emsr::detail::gamma(a);

      return U1 + U2;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    tricomi_u_asymp(_Tp a, _Tp c, _Tp z)
    {
      const auto s_eps = emsr::epsilon(z);
      const unsigned int s_max_iter = 100000u;
      auto b = a - c + _Tp{1};
      auto term = _Tp{1};
      auto _Usum = _Tp{1};
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  term *= -(a + _Tp(k - 1)) * (b + _Tp(k - 1))
		  / _Tp(k) / z;
	  _Usum += term;
	  if (std::abs(term) < s_eps * std::abs(_Usum))
	    break;
	}
      return std::pow(z, -a) * _Usum;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    tricomi_u_c_pos_int(_Tp a, int m, _Tp z)
    {
      const auto s_eps = emsr::epsilon(z);
      const unsigned int s_max_iter = 100000u;
//std::cout << '\n';

      auto term1 = _Tp{1};
      auto _U1 = term1;
      for (auto k = 1; k <= m - 2; ++k)
	{
	  term1 *= (a + _Tp(-m + k)) * z
		   / _Tp(1 - m + k) / _Tp(k);
	  _U1 += term1;
//std::cout << "_U1 = " << _U1 << '\n';
	}
      _U1 *= emsr::detail::factorial<_Tp>(m - 2)
	   / emsr::detail::gamma(a)
	   / std::pow(z, _Tp(m - 1));

      const auto b = a + _Tp(1 - m);
      auto psi2 = emsr::detail::digamma(a)
		  - emsr::detail::digamma<_Tp>(1)
		  - emsr::detail::digamma<_Tp>(m)
		  + std::log(z);
      auto fact2 = _Tp{1};
      auto _U2 = psi2;
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  psi2 += _Tp{1} / (a + _Tp(k - 1))
		  - _Tp{1} / _Tp(k)
		  - _Tp{1} / _Tp(m + k - 1);
	  fact2 *= (a + _Tp(k)) * z / _Tp(m + k) / _Tp(k);
	  auto term2 = psi2 * fact2;
	  _U2 += term2;
//std::cout << "_U2 = " << _U2 << '\n';
	  if (std::abs(term2) < s_eps * std::abs(_U2))
	    break;
	}
      _U2 *= ((m & 1) ? -1 : +1)
	   / emsr::detail::factorial<_Tp>(m - 1)
	   / emsr::detail::gamma(b);

      return _U1 + _U2;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    tricomi_u_c_nonpos_int(_Tp a, _Tp m, _Tp z)
    {
      auto b = a + _Tp(1 - m);
      return std::pow(z, _Tp(1 - m))
	   * tricomi_u_c_pos_int(b, 2 - m, z);
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    tricomi_u_ac_int(_Tp n, _Tp m, _Tp z)
    {
      const unsigned int s_max_iter = 100000u;
      const auto s_eps = emsr::epsilon(z);
      auto term = _Tp{1};
      auto _Usum = _Tp{1};
      for (auto k = 1u; k < -n; ++k)
	{
	  term *= _Tp(n + k - 1) / _Tp(m + k - 1) / _Tp(k) * z;
	  _Usum += term;
	  if (std::abs(term) < s_eps * std::abs(_Usum))
	    break;
	}
      return ((n & 1) ? -1 : +1)
	   * emsr::rising_factorial(m, -n) * _Usum;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    tricomi_u_series(_Tp a, _Tp c, _Tp z)
    {
      return 0;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    tricomi_u(_Tp a, _Tp c, _Tp z)
    {
      auto aint = emsr::fp_is_integer(a);
      auto cint = emsr::fp_is_integer(c);
      if (cint)
	{
	  if (cint() > 0)
	    return tricomi_u_c_pos_int(a, cint(), z);
	  else
	    return tricomi_u_c_nonpos_int(a, cint(), z);
	}
      else
	return tricomi_u_series(a, c, z);
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    whittaker_m(_Tp kappa, _Tp mu, _Tp z)
    {
      return std::exp(-z / 2) * std::pow(z, 0.5 + mu)
	   * emsr::conf_hyperg(0.5 + mu - kappa, 1 + 2 * mu, z);
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    whittaker_w(_Tp kappa, _Tp mu, _Tp z)
    {
      return std::exp(-z / 2) * std::pow(z, 0.5 + mu)
	   * tricomi_u(0.5 + mu - kappa, 1 + 2 * mu, z);
    }

//} // namespace detail
//} // namespace emsr

template<typename _Tp>
  void
  test_tricomi_u(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto a = _Tp{6} / _Tp{5};
    const auto c = _Tp{1} / _Tp{5};
    const auto del = _Tp{1} / _Tp{10};
    for (int i = 0; i < +200; ++i)
    {
      auto z = del * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << tricomi_u_naive(a, c, z)
		<< '\n';
    }

    std::cout << "\nInteger c = m\n";
    const auto z = _Tp{1} / _Tp{2};
    std::cout << " a = " << std::setw(6) << a << '\n';
    std::cout << " z = " << std::setw(6) << z << '\n';
    for (auto m = 1u; m <= +20; ++m)
    {
      std::cout << ' ' << std::setw(6) << m
		<< ' ' << std::setw(width) << tricomi_u_naive(a, _Tp(m), z)
		<< ' ' << std::setw(width) << tricomi_u_c_pos_int(a, m, z)
		<< '\n';
    }
  }

int
main()
{
  test_tricomi_u(1.0);
}
