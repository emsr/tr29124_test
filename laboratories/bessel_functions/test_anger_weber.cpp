/**
 *
 */

#include <cmath>
#include <emsr/math_util.h>
#include <emsr/numeric_limits.h>
#include <emsr/specfun.h>

  /**
   * 
   */
  template<typename _Tp>
    struct anger_weber_t
    {
      _Tp nu;
      _Tp z;
      _Tp J_value;
      _Tp E_value;
    };

  /**
   * A smart Gamma reciprocal function/iterator.
   */
  template<typename _Tp>
    struct _GammaReciprocal
    {
      explicit _GammaReciprocal(_Tp a)
      : _M_arg(a),
        _M_int(emsr::fp_is_integer(a))
      { }

      _Tp
      operator()()
      {
	if (this->_M_int && this->_M_arg <= _Tp{0})
	  return _Tp{0};
	else
	  {
	    if (this->_M_start)
	      return (this->_M_gam /= (this->_M_arg += _Tp{1}));
	    else
	      {
	        this->_M_start = true;
	        this->_M_gam = emsr::detail::gamma_reciprocal(this->_M_arg);
		return this->_M_gam;
	      }
	  }
      }

      _GammaReciprocal&
      operator++()
      {
	this->_M_arg += _Tp{1};
	return *this;
      }

      _GammaReciprocal&
      operator++(int)
      {
	auto temp = *this;
	this->_M_arg += _Tp{1};
	return temp;
      }

      _Tp _M_arg;
      bool _M_int;
      bool _M_start = false;
    };

  /**
   * Compute a sum used for Anger and Weber functions:
   * @f[
   *    S_1(\nu,z) = \sum_{k=0}^{\infty}\frac{(-1)^k(\frac{z}{2})^{2k}}
   *   {\Gamma(k+\frac{\nu}{2}+1) \Gamma(k-\frac{\nu}{2}+1)}
   * @f]
   */
  template<typename _Tp>
    _Tp
    anger_weber_sum_1(_Tp nu, _Tp z)
    {
      const auto s_max_iter = 10000u;
      const auto s_eps = emsr::epsilon(z);
      const auto z2 = z / _Tp{2};
      auto _GamArg11 = _Tp{1} + nu / _Tp{2};
      auto _GamArg12 = _Tp{1} - nu / _Tp{2};
      auto _Gam11 = std::tgamma(_GamArg11);
      auto _Gam12 = std::tgamma(_GamArg12);
      auto term1 = _Tp{1} / (_Gam11 * _Gam12);
      auto _S1 = term1;
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  term1 *= -z2 / _GamArg11 * z2 / _GamArg12;
	  _S1 += term1;
	  _GamArg11 += _Tp{1};
	  _GamArg12 += _Tp{1};

	  if (std::abs(term1) < s_eps * std::abs(_S1))
	    return _S1;
	}
      return _Tp{0};
    }

  /**
   * Compute a sum used for Anger and Weber functions.
   * Assumes n == 2m, m > 0.
   * @f[
   *    S_1(\nu,z) = \sum_{k=0}^{\infty}\frac{(-1)^k(\frac{z}{2})^{2k}}
   *   {\Gamma(k+\frac{\nu}{2}+1) \Gamma(k-\frac{\nu}{2}+1)}
   * @f]
   */
  template<typename _Tp>
    _Tp
    anger_weber_sum_1_even_int(int n, _Tp z)
    {
      const auto s_max_iter = 10000u;
      const auto s_eps = emsr::epsilon(z);
      const auto z2 = z / _Tp{2};
      const auto m = n / 2;
      const auto k_start = m;
      auto _GamArg11 = _Tp(1 + m + k_start);
      auto _Gam11 = std::tgamma(_GamArg11);
      auto _GamArg12 = _Tp(1 - m + k_start);
      auto _Gam12 = std::tgamma(_GamArg12);
      auto term1 = ((k_start & 1) ? -1 : +1)
		   * std::pow(z2, _Tp(2 * k_start)) / (_Gam11 * _Gam12);
      auto _S1 = term1;
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  term1 *= -z2 / _GamArg11 * z2 / _GamArg12;
	  _S1 += term1;
	  _GamArg11 += _Tp{1};
	  _GamArg12 += _Tp{1};

	  if (std::abs(term1) < s_eps * std::abs(_S1))
	    return _S1;
	}
      return _Tp{0};
    }

  /**
   * Compute a sum used for Anger and Weber functions:
   * @f[
   *    S_2(\nu,z) = \sum_{k=0}^{\infty}\frac{(-1)^k(\frac{z}{2})^{2k+1}}
   *   {\Gamma(k+\frac{\nu}{2}+\frac{3}{2}) \Gamma(k-\frac{\nu}{2}+\frac{3}{2})}
   * @f]
   */
  template<typename _Tp>
    _Tp
    anger_weber_sum_2(_Tp nu, _Tp z)
    {
      const auto s_max_iter = 10000u;
      const auto s_eps = emsr::epsilon(z);
      const auto z2 = z / _Tp{2};
      auto _GamArg21 = _Tp{3} / _Tp{2} + nu / _Tp{2};
      auto _GamArg22 = _Tp{3} / _Tp{2} - nu / _Tp{2};
      auto _Gam21 = std::tgamma(_GamArg21);
      auto _Gam22 = std::tgamma(_GamArg22);
      auto term2 = z2 / (_Gam21 * _Gam22);
      auto _S2 = term2;
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  term2 *= -z2 / _GamArg21 * z2 / _GamArg22;
	  _S2 += term2;
	  _GamArg21 += _Tp{1};
	  _GamArg22 += _Tp{1};

	  if (std::abs(term2) < s_eps * std::abs(_S2))
	    return _S2;
	}
      return _Tp{0};
    }

  /**
   * Compute a sum used for Anger and Weber functions.
   * Assumes n == 2m+1, m > 0.
   * @f[
   *    S_2(\nu,z) = \sum_{k=0}^{\infty}\frac{(-1)^k(\frac{z}{2})^{2k+1}}
   *   {\Gamma(k+\frac{\nu}{2}+\frac{3}{2}) \Gamma(k-\frac{\nu}{2}+\frac{3}{2})}
   * @f]
   */
  template<typename _Tp>
    _Tp
    anger_weber_sum_2_odd_int(int n, _Tp z)
    {
      const auto s_max_iter = 10000u;
      const auto s_eps = emsr::epsilon(z);
      const auto z2 = z / _Tp{2};
      const auto m = (n - 1) / 2;
      const auto k_start = m;
      auto _GamArg21 = _Tp{2} + m + k_start;
      auto _Gam21 = std::tgamma(_GamArg21);
      auto _GamArg22 = _Tp{1} - m + k_start;
      auto _Gam22 = std::tgamma(_GamArg22);
      auto term2 = ((k_start & 1) ? -1 : +1)
		   * std::pow(z2, _Tp(2 * k_start + 1)) / (_Gam21 * _Gam22);
      auto _S2 = term2;
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  term2 *= -z2 / _GamArg21 * z2 / _GamArg22;
	  _S2 += term2;
	  _GamArg21 += _Tp{1};
	  _GamArg22 += _Tp{1};

	  if (std::abs(term2) < s_eps * std::abs(_S2))
	    return _S2;
	}
      return _Tp{0};
    }

  /**
   * 
   */
  template<typename _Tp>
    anger_weber_t<_Tp>
    anger_weber_sum_new(_Tp nu, _Tp z)
    {
      if (nu < _Tp{0})
	{
	  auto AW = anger_weber_sum_new(-nu, -z);
	  AW.E_value = -AW.E_value;
	  return AW;
	}
      else
	{
	  auto nuint = emsr::fp_is_integer(nu);

	  auto _S1 = _Tp{0};
	  if (nuint && nuint() > 0 && nuint() % 2 == 0)
	    _S1 = anger_weber_sum_1_even_int(nuint(), z);
	  else
	    _S1 = anger_weber_sum_1(nu, z);

	  auto _S2 = _Tp{0};
	  if (nuint && nuint() > 0 && nuint() % 2 == 1)
	    _S2 = anger_weber_sum_2_odd_int(nuint(), z);
	  else
	    _S2 = anger_weber_sum_2(nu, z);

	  auto ph = emsr::detail::sincos_pi(nu / _Tp{2});
	  return anger_weber_t<_Tp>{nu, z,
				      ph.cos_v * _S1
				    + ph.sin_v * _S2,
				      ph.sin_v * _S1
				    - ph.cos_v * _S2};
	}
    }

  /**
   * 
   */
  template<typename _Tp>
    anger_weber_t<_Tp>
    anger_weber_sum(_Tp nu, _Tp z)
    {
      const auto s_eps = emsr::epsilon(z);

      auto nuint = emsr::fp_is_integer(nu);

      if (nu < _Tp{0})
	{
	  auto AW = anger_weber_sum(-nu, -z);
	  AW.E_value = -AW.E_value;
	  return AW;
	}
      else if (nuint && nuint() > 1)
	{
	  auto n = nuint();
	  if (n & 1)
	    {
	      const auto z2 = z / _Tp{2};
	      auto _GamArg11 = _Tp{1} + nu / _Tp{2};
	      auto _GamArg12 = _Tp{1} - nu / _Tp{2};
	      auto _Gam11 = std::tgamma(_GamArg11);
	      auto _Gam12 = std::tgamma(_GamArg12);
	      auto term1 = _Tp{1} / (_Gam11 * _Gam12);
	      auto _S1 = term1;
	      for (auto k = 1u; k < 10000u; ++k)
		{
		  term1 *= -z2 / _GamArg11 * z2 / _GamArg12;
		  _S1 += term1;
		  _GamArg11 += _Tp{1};
		  _GamArg12 += _Tp{1};

		  if (std::abs(term1) < s_eps * std::abs(_S1))
		    break;
		}
	      return anger_weber_t<_Tp>{nu, z,
					  _Tp{0},
					  -(((n / 2) & 1) ? -1 : +1) * _S1};
	    }
	  else
	    {
	      const auto z2 = z / _Tp{2};
	      auto _GamArg21 = _Tp{3} / _Tp{2} + nu / _Tp{2};
	      auto _GamArg22 = _Tp{3} / _Tp{2} - nu / _Tp{2};
	      auto _Gam21 = std::tgamma(_GamArg21);
	      auto _Gam22 = std::tgamma(_GamArg22);
	      auto term2 = z2 / (_Gam21 * _Gam22);
	      auto _S2 = term2;
	      for (auto k = 1u; k < 10000u; ++k)
		{
		  term2 *= -z2 / _GamArg21 * z2 / _GamArg22;
		  _S2 += term2;
		  _GamArg21 += _Tp{1};
		  _GamArg22 += _Tp{1};

		  if (std::abs(term2) < s_eps * std::abs(_S2))
		    break;
		}
	      return anger_weber_t<_Tp>{nu, z,
					  _Tp{0},
					 -((n / 2) & 1 ? -1 : +1) * _S2};
	    }
	}
      else
	{
	  const auto z2 = z / _Tp{2};
	  auto _GamArg11 = _Tp{1} + nu / _Tp{2};
	  auto _GamArg12 = _Tp{1} - nu / _Tp{2};
	  auto _GamArg21 = _Tp{3} / _Tp{2} + nu / _Tp{2};
	  auto _GamArg22 = _Tp{3} / _Tp{2} - nu / _Tp{2};
	  auto _Gam11 = std::tgamma(_GamArg11);
	  auto _Gam12 = std::tgamma(_GamArg12);
	  auto _Gam21 = std::tgamma(_GamArg21);
	  auto _Gam22 = std::tgamma(_GamArg22);
	  auto term1 = _Tp{1} / (_Gam11 * _Gam12);
	  auto _S1 = term1;
	  auto term2 = z2 / (_Gam21 * _Gam22);
	  auto _S2 = term2;
	  for (auto k = 1u; k < 10000u; ++k)
	    {
	      term1 *= -z2 / _GamArg11 * z2 / _GamArg12;
	      _S1 += term1;
	      _GamArg11 += _Tp{1};
	      _GamArg12 += _Tp{1};

	      term2 *= -z2 / _GamArg21 * z2 / _GamArg22;
	      _S2 += term2;
	      _GamArg21 += _Tp{1};
	      _GamArg22 += _Tp{1};

	      if (std::abs(term1) < s_eps * std::abs(_S1)
	       && std::abs(term2) < s_eps * std::abs(_S2))
		break;
	    }
	  //auto [sin, cos] = sincos_pi(nu / _Tp{2});
	  auto ph = emsr::detail::sincos_pi(nu / _Tp{2});
	  return anger_weber_t<_Tp>{nu, z,
				      ph.cos_v * _S1
				    + ph.sin_v * _S2,
				      ph.sin_v * _S1
				    - ph.cos_v * _S2};
	}
    }

  /**
   * Compute Anger and Weber functions for fixed order @f$ \nu @f$
   * and large agument @f$ |z| @f$.
   *
   * @see http://dlmf.nist.gov/11.10#i
   */
  template<typename _Tp>
    anger_weber_t<_Tp>
    anger_weber_asymp_arg(_Tp nu, _Tp z)
    {
      using _Real = decltype(std::real(z));
      const auto s_eps = emsr::epsilon<_Real>();
      const auto s_pi = emsr::pi_v<_Real>;
      const auto s_max_iter = 1000u;
      const auto z2 = z * z;

      auto F_z2k = _Tp{1};
      auto G_z2k = _Tp{1};

      auto Fsum = F_z2k;
      auto Gsum = G_z2k;
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  F_z2k *= (nu - _Tp(2 * k - 1)) * (nu + _Tp(2 * k - 1))
		   / z2;
	  Fsum += F_z2k;
	  G_z2k *= (nu - _Tp(2 * k)) * (nu + _Tp(2 * k))
		   / z2;
	  Gsum += G_z2k;
	}

      auto ph = emsr::detail::sincos_pi(nu / _Tp{2});
      auto _Bess = cyl_bessel(nu, z);
      return anger_weber_t<_Tp>{nu, z,
				  _Bess._J_value
				    + ph.sin_v
				* (Fsum + nu * Gsum / z) / s_pi / z,
				 -_Bess._N_value
				    - (_Tp{1} + ph.cos_v) * Fsum
					/ s_pi / z
				    - (_Tp{1} - ph.cos_v) * Gsum
					* nu / s_pi / z / z};
    }

  /**
   * Compute Anger and Weber functions for large order @f$ |\nu| @f$
   * and fixed agument @f$ z @f$.
   *
   * @see http://dlmf.nist.gov/11.10#ii
   */
  template<typename _Tp>
    anger_weber_t<_Tp>
    anger_weber_asymp_order(_Tp nu, _Tp z)
    {
      using _Real = decltype(std::real(z));
      const auto s_pi = emsr::pi_v<_Real>;
      const auto sinnp = emsr::sin_pi(nu);
      const auto sinnpd2 = emsr::sin_pi(nu / _Tp{2});
      const auto cosnpd2 = emsr::cos_pi(nu / _Tp{2});
      const auto nufact = nu * z / (nu * nu - _Tp{1});
      return anger_weber_t<_Tp>{nu, z,
				  sinnp * (_Tp{1} - nufact) / nu / s_pi,
				  _Tp{2} * (sinnpd2 + nufact * cosnpd2)
					 / nu / s_pi};
    }

  /**
   * Compute Anger and Weber functions for large order @f$ \nu @f$
   * and fixed ratio @f$ z/\nu @f$.
   *
   * @see http://dlmf.nist.gov/11.10#iii
   */
  template<typename _Tp>
    anger_weber_t<_Tp>
    anger_weber_asymp_uniform(_Tp nu, _Tp z)
    {
    }

  /**
   * Use the reciprocal gamma function... Fails.  WTF.
   */
  template<typename _Tp>
    anger_weber_t<_Tp>
    anger_weber_sum_recip(_Tp nu, _Tp z)
    {
      //using _Val = _Tp;
      //using _Real = emsr::num_traits_t<_Val>;
      const auto s_eps = emsr::epsilon(z);

      const auto z2 = z / _Tp{2};
      auto _GamArg11 = _Tp{1} + nu / _Tp{2};
      auto _GamArg12 = _Tp{1} - nu / _Tp{2};
      auto _GamArg21 = _Tp{3} / _Tp{2} + nu / _Tp{2};
      auto _GamArg22 = _Tp{3} / _Tp{2} - nu / _Tp{2};
      auto _Gam11 = emsr::detail::gamma_reciprocal(_GamArg11);
      auto _Gam12 = emsr::detail::gamma_reciprocal(_GamArg12);
      auto _Gam21 = emsr::detail::gamma_reciprocal(_GamArg21);
      auto _Gam22 = emsr::detail::gamma_reciprocal(_GamArg22);
      auto term1 = _Gam11 * _Gam12;
      auto _S1 = term1;
      auto term2 = z2 * _Gam21 * _Gam22;
      auto _S2 = term2;
      for (auto k = 1u; k < 10000u; ++k)
	{
	  term1 *= -z2 / _GamArg11 * z2 / _GamArg12;
	  _S1 += term1;
	  _GamArg11 += _Tp{1};
	  _GamArg12 += _Tp{1};

	  term2 *= -z2 / _GamArg21 * z2 / _GamArg22;
	  _S2 += term2;
	  _GamArg21 += _Tp{1};
	  _GamArg22 += _Tp{1};

	  if (std::abs(term1) < s_eps * std::abs(_S1)
	   && std::abs(term2) < s_eps * std::abs(_S2))
	    break;
	}
      //auto [sin, cos] = sincos_pi(nu / _Tp{2});
      auto ph = emsr::detail::sincos_pi(nu / _Tp{2});
      return anger_weber_t<_Tp>{nu, z,
				  ph.cos_v * _S1
				+ ph.sin_v * _S2,
				  ph.sin_v * _S1
				- ph.cos_v * _S2};
    }

  /**
   * Compute the Anger @f$ {\boldmath J}_\nu(z) @f$
   * and Weber @f$ {\boldmath E}_\nu(z) @f$ functions
   * for order @f$ \nu @f$ and agument @f$ z @f$.
   *
   * @see http://dlmf.nist.gov/11.10#ii
   */
  template<typename _Tp>
    _Tp
    assoc_anger_weber_asymp(_Tp nu, _Tp z)
    {
      auto _Bessel = cyl_bessel(nu, z);
      auto _Weber = anger_weber(nu, z);
    }

  /**
   * Compute the associated Anger and Weber function @f$ A_\nu(z) @f$
   * for order @f$ \nu @f$ and agument @f$ z @f$.
   * We use the relationship:
   * @f[
   *   \boldmath{J}_\nu(z) = J_\nu(z) + \sin(\nu \pi) \boldmath{A}_\nu(z)
   * @f]
   * Note that for integer order n @f$ \boldmath{A}_n(z) = 0 @f$
   * since @f$ \boldmath{J}_n(z) = J_n(z) @f$.
   *
   * @see http://dlmf.nist.gov/11.10#v
   */
  template<typename _Tp>
    _Tp
    assoc_anger_weber(_Tp nu, _Tp z)
    {
      auto _Bessel = emsr::cyl_bessel_j(nu, z);
      auto _Weber = anger_weber_sum_new(nu, z);
      if (emsr::fp_is_integer(nu))
	return _Tp{0};
      else
        return (_Weber.J_value - _Bessel) / emsr::sin_pi(nu);
    }

  /**
   * @todo: Find out if
   * @f[
   *   \boldmath{J}_\nu(z) + i\boldmath{E}_\nu(z)
   * @f]
   * is a thing.
   */

  /**
   * \frac{1 - \cos(\pi x)}{\pi x}
   */
  template<typename _Tp>
  _Tp
  cosc_pi(_Tp x)
  {
    const auto s_eps = emsr::epsilon(x);
    const auto s_pi = emsr::pi_v<_Tp>;
    if (std::abs(x) < _Tp{100} * s_eps)
      return s_pi * x / _Tp{2};
    else
      return (_Tp{1} - emsr::cos_pi(x)) / s_pi / x;
  }


template<typename _Tp>
  void
  test_anger_weber(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    //std::cout << "\n\n Write J and E values\n";
    //std::cout << " --------------------\n";
    const auto twk = _Tp{1}/_Tp{1000};
    for (auto nu : {_Tp{0}, _Tp{1}/_Tp{2}, _Tp{1}, _Tp{3}/_Tp{2},
		    _Tp{2} - twk, _Tp{2},
		    _Tp{3} - twk, _Tp{3},
		    _Tp{5}})
      {
	std::cout << "\n\n nu = " << std::setw(4) << nu << '\n';
	std::cout << ' ' << std::setw(4) << "z"
		  << ' ' << std::setw(width) << "Jbold"
		  << ' ' << std::setw(width) << "Ebold"
		  << '\n';
	std::cout << ' ' << std::setw(4) << "-"
		  << ' ' << std::setw(width) << "-----"
		  << ' ' << std::setw(width) << "-----"
		  << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (int k = -80; k <= 80; ++k)
	  {
	    auto z = del * k;
	    //auto AW = anger_weber_sum(nu, z);
	    auto AW = anger_weber_sum_new(nu, z);
	    std::cout << ' ' << std::setw(4) << AW.z
		      << ' ' << std::setw(width) << AW.J_value
		      << ' ' << std::setw(width) << AW.E_value
		      << '\n';
	  }
      }

    std::cout << "\n\n Test J and E values at zero\n";
    std::cout << " ---------------------------\n";
    for (auto nu : {_Tp{0}, _Tp{1}/_Tp{2}, _Tp{1}, _Tp{3}/_Tp{2},
		    _Tp{2} - twk, _Tp{2},
		    _Tp{3} - twk, _Tp{3},
		    _Tp{5}})
      {
	auto AW = anger_weber_sum_new(nu, _Tp{0});
	std::cout << "\n\n nu = " << std::setw(4) << AW.nu << '\n';
	std::cout << ' ' << std::setw(width) << AW.J_value
		  << ' ' << std::setw(width) << emsr::sinc_pi(nu)
		  << ' ' << std::setw(width) << AW.E_value
		  << ' ' << std::setw(width) << cosc_pi(nu)
		  << '\n';
      }

    std::cout << "\n\n Test J values for integer order\n";
    std::cout << " -------------------------------\n";
    for (auto nu : {_Tp{0}, _Tp{1}, _Tp{2}, _Tp{3}, _Tp{5}})
      {
	std::cout << "\n\n nu = " << std::setw(4) << nu << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (int k = 0; k <= 80; ++k)
	  {
	    auto z = del * k;
	    auto AW = anger_weber_sum_new(nu, z);
	    std::cout << ' ' << std::setw(4) << AW.z
		      << ' ' << std::setw(width) << AW.J_value
		      << ' ' << std::setw(width) << emsr::cyl_bessel_j(nu, z)
		      << '\n';
	  }
      }

    std::cout << "\n\n Write A values\n";
    std::cout << " -------------------------------\n";
    for (auto nu : {_Tp{0}, _Tp{1}/_Tp{2}, _Tp{1}, _Tp{3}/_Tp{2},
		    _Tp{2} - twk, _Tp{2},
		    _Tp{3} - twk, _Tp{3},
		    _Tp{5}})
      {
	std::cout << "\n\n nu = " << std::setw(4) << nu << '\n';
	std::cout << ' ' << std::setw(4) << "z"
		  << ' ' << std::setw(width) << "Abold"
		  << '\n';
	std::cout << ' ' << std::setw(4) << "-"
		  << ' ' << std::setw(width) << "-----"
		  << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (int k = 0; k <= 80; ++k)
	  {
	    auto z = del * k;
	    auto AW = assoc_anger_weber(nu, z);
	    std::cout << ' ' << std::setw(4) << z
		      << ' ' << std::setw(width) << AW
		      << '\n';
	  }
      }

    std::cout << "\n\n Test A values for integer order\n";
    std::cout << " -------------------------------\n";
    for (auto nu : {_Tp{0}, _Tp{1}, _Tp{2}, _Tp{3}, _Tp{5}})
      {
	std::cout << "\n\n nu = " << std::setw(4) << nu << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (int k = 0; k <= 80; ++k)
	  {
	    auto z = del * k;
	    auto AWm = assoc_anger_weber(nu - twk, z);
	    auto AW = assoc_anger_weber(nu, z);
	    auto AWp = assoc_anger_weber(nu + twk, z);
	    std::cout << ' ' << std::setw(4) << z
		      << ' ' << std::setw(width) << AWm
		      << ' ' << std::setw(width) << AW
		      << ' ' << std::setw(width) << AWp
		      << '\n';
	  }
      }
  }

int
main()
{
  auto AW2 [[maybe_unused]] = anger_weber_sum_new(2.0, -8.0);
  auto BX2 [[maybe_unused]] = anger_weber_sum_new(1.999, -8.0);
  auto AW3 [[maybe_unused]] = anger_weber_sum_new(3.0, -8.0);
  auto BX3 [[maybe_unused]] = anger_weber_sum_new(2.999, -8.0);
  test_anger_weber(1.0);
}
