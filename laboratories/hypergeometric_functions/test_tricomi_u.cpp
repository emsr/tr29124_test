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
  template<typename Tp>
    Tp
    tricomi_u_naive(Tp a, Tp c, Tp x)
    {
      auto U1 = Tp{};
      auto b = a - c + Tp{1};
      auto ib = emsr::fp_is_integer(b);
      if (!ib || ib() > 0)
	U1 = emsr::detail::gamma(Tp{1} - c)
	     * emsr::detail::conf_hyperg(a, c, x)
	     / emsr::detail::gamma(b);

      auto U2 = Tp{};
      auto ia = emsr::fp_is_integer(a);
      if (!ia || ia() > 0)
	U2 = emsr::detail::gamma(c - Tp{1})
	     * std::pow(x, Tp{1} - c)
	     * emsr::detail::conf_hyperg(b, Tp{2} - c, x)
	     / emsr::detail::gamma(a);

      return U1 + U2;
    }

  /**
   * 
   */
  template<typename Tp>
    Tp
    tricomi_u_asymp(Tp a, Tp c, Tp z)
    {
      const auto s_eps = emsr::epsilon(z);
      const unsigned int s_max_iter = 100000u;
      auto b = a - c + Tp{1};
      auto term = Tp{1};
      auto _Usum = Tp{1};
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  term *= -(a + Tp(k - 1)) * (b + Tp(k - 1))
		  / Tp(k) / z;
	  _Usum += term;
	  if (std::abs(term) < s_eps * std::abs(_Usum))
	    break;
	}
      return std::pow(z, -a) * _Usum;
    }

  /**
   * 
   */
  template<typename Tp>
    Tp
    tricomi_u_c_pos_int(Tp a, int m, Tp z)
    {
      const auto s_eps = emsr::epsilon(z);
      const unsigned int s_max_iter = 100000u;
//std::cout << '\n';

      auto term1 = Tp{1};
      auto _U1 = term1;
      for (auto k = 1; k <= m - 2; ++k)
	{
	  term1 *= (a + Tp(-m + k)) * z
		   / Tp(1 - m + k) / Tp(k);
	  _U1 += term1;
//std::cout << "_U1 = " << _U1 << '\n';
	}
      _U1 *= emsr::detail::factorial<Tp>(m - 2)
	   / emsr::detail::gamma(a)
	   / std::pow(z, Tp(m - 1));

      const auto b = a + Tp(1 - m);
      auto psi2 = emsr::detail::digamma(a)
		  - emsr::detail::digamma<Tp>(1)
		  - emsr::detail::digamma<Tp>(m)
		  + std::log(z);
      auto fact2 = Tp{1};
      auto _U2 = psi2;
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  psi2 += Tp{1} / (a + Tp(k - 1))
		  - Tp{1} / Tp(k)
		  - Tp{1} / Tp(m + k - 1);
	  fact2 *= (a + Tp(k)) * z / Tp(m + k) / Tp(k);
	  auto term2 = psi2 * fact2;
	  _U2 += term2;
//std::cout << "_U2 = " << _U2 << '\n';
	  if (std::abs(term2) < s_eps * std::abs(_U2))
	    break;
	}
      _U2 *= ((m & 1) ? -1 : +1)
	   / emsr::detail::factorial<Tp>(m - 1)
	   / emsr::detail::gamma(b);

      return _U1 + _U2;
    }

  /**
   * 
   */
  template<typename Tp>
    Tp
    tricomi_u_c_nonpos_int(Tp a, Tp m, Tp z)
    {
      auto b = a + Tp(1 - m);
      return std::pow(z, Tp(1 - m))
	   * tricomi_u_c_pos_int(b, 2 - m, z);
    }

  /**
   * 
   */
  template<typename Tp>
    Tp
    tricomi_u_ac_int(Tp n, Tp m, Tp z)
    {
      const unsigned int s_max_iter = 100000u;
      const auto s_eps = emsr::epsilon(z);
      auto term = Tp{1};
      auto _Usum = Tp{1};
      for (auto k = 1u; k < -n; ++k)
	{
	  term *= Tp(n + k - 1) / Tp(m + k - 1) / Tp(k) * z;
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
  template<typename Tp>
    Tp
    tricomi_u_series(Tp a, Tp c, Tp z)
    {
      return 0;
    }

  /**
   * 
   */
  template<typename Tp>
    Tp
    tricomi_u(Tp a, Tp c, Tp z)
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
  template<typename Tp>
    Tp
    whittaker_m(Tp kappa, Tp mu, Tp z)
    {
      return std::exp(-z / 2) * std::pow(z, 0.5 + mu)
	   * emsr::conf_hyperg(0.5 + mu - kappa, 1 + 2 * mu, z);
    }

  /**
   * 
   */
  template<typename Tp>
    Tp
    whittaker_w(Tp kappa, Tp mu, Tp z)
    {
      return std::exp(-z / 2) * std::pow(z, 0.5 + mu)
	   * tricomi_u(0.5 + mu - kappa, 1 + 2 * mu, z);
    }

//} // namespace detail
//} // namespace emsr

template<typename Tp>
  void
  test_tricomi_u(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto a = Tp{6} / Tp{5};
    const auto c = Tp{1} / Tp{5};
    const auto del = Tp{1} / Tp{10};
    for (int i = 0; i < +200; ++i)
    {
      auto z = del * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << tricomi_u_naive(a, c, z)
		<< '\n';
    }

    std::cout << "\nInteger c = m\n";
    const auto z = Tp{1} / Tp{2};
    std::cout << " a = " << std::setw(6) << a << '\n';
    std::cout << " z = " << std::setw(6) << z << '\n';
    for (auto m = 1u; m <= +20; ++m)
    {
      std::cout << ' ' << std::setw(6) << m
		<< ' ' << std::setw(width) << tricomi_u_naive(a, Tp(m), z)
		<< ' ' << std::setw(width) << tricomi_u_c_pos_int(a, m, z)
		<< '\n';
    }
  }

int
main()
{
  test_tricomi_u(1.0);
}
