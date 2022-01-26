/**
 *
 */

/**
 * Look at the formula for the reciprocal of the gamma for the Temme gamma
 * \frac{1}{\Gamma(1 +- \mu)}
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

#include <emsr/float128_io.h>

//#include <mpreal.h>
//#include <emsr/math_const_mpreal.h>
//#include <math_mpreal.h>
//#include <emsr/numeric_limits_mpreal.h>

  /**
   * 
   */
  template<typename _Tp>
    std::vector<emsr::num_traits_t<_Tp>>
    gamma_reciprocal_series_coef(std::size_t n, _Tp proto = _Tp{})
    {
      using _Val = _Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto s_eps = emsr::epsilon(std::real(proto));
      const auto s_gamma_e = emsr::egamma_v<_Tp>;
      auto sign = [](std::size_t i){ return (i & 1u) == 1u ? -1 : +1; };
      std::vector<_Real> c;
      c.push_back(_Real{0});
      c.push_back(_Real{1});
      for (auto j = 1u; j < n; ++j)
	{
	  auto sum = _Real{0};
	  for (auto k = 1u; k < j; ++k)
	    sum += sign(k) * c[k]
		   * (_Real{1} + emsr::detail::riemann_zeta_m_1(_Real(j + 1 - k)));
	  c.push_back((s_gamma_e * c[j] + sign(j) * sum) / j);
	  if (std::abs(c.back()) < s_eps)
	    break;
	}
      return c;
    }

  /**
   * Return the reciprocal of the Gamma function by series.
   * The reciprocal of the Gamma function is given by
   * @f[
   *   \frac{1}{\Gamma(a)} = \sum_{k=1}^{\infty} c_k a^k
   * @f]
   * where the coefficients are defined by recursion:
   * @f[
   *   c_{k+1} = \frac{1}{k}\left[\gamma_E c_k
   *           + (-1)^k\sum_{j=1}^{k-1}(-1)^j\zeta(j+1-k)c_j\right]
   * @f]
   * where @f$ c_1 = 1 @f$
   */
  template<typename _Tp>
    _Tp
    gamma_reciprocal_series(_Tp a)
    {
      static constexpr std::array<__float128, 50>
      s_c
      {{
	 0.0000000000000000000000000000000000000000Q,
	 1.0000000000000000000000000000000000000000Q,
	 0.5772156649015328606065120900824024310432Q,
	-0.6558780715202538810770195151453904812811Q,
	-0.0420026350340952355290039348754298187119Q,
	 0.1665386113822914895017007951021052357187Q,
	-0.0421977345555443367482083012891873913015Q,
	-0.0096219715278769735621149216723481989747Q,
	 0.0072189432466630995423950103404465727093Q,
	-0.0011651675918590651121139710840183886674Q,
	-0.0002152416741149509728157299630536478063Q,
	 0.0001280502823881161861531986263281643238Q,
	-0.0000201348547807882386556893914210218186Q,
	-0.0000012504934821426706573453594738330926Q,
	 0.0000011330272319816958823741296203307448Q,
	-0.0000002056338416977607103450154130020573Q,
	 0.0000000061160951044814158178624986828556Q,
	 0.0000000050020076444692229300556650480601Q,
	-0.0000000011812745704870201445881265654365Q,
	 0.0000000001043426711691100510491540332313Q,
	 0.0000000000077822634399050712540499373115Q,
	-0.0000000000036968056186422057081878158781Q,
	 0.0000000000005100370287454475979015481319Q,
	-0.0000000000000205832605356650678322242954Q,
	-0.0000000000000053481225394230179823700171Q,
	 0.0000000000000012267786282382607901588941Q,
	-0.0000000000000001181259301697458769513765Q,
	 0.0000000000000000011866922547516003325796Q,
	 0.0000000000000000014123806553180317815559Q,
	-0.0000000000000000002298745684435370206591Q,
	 0.0000000000000000000171440632192733743337Q,
	 0.0000000000000000000001337351730493693114Q,
	-0.0000000000000000000002054233551766672789Q,
	 0.0000000000000000000000273603004860799984Q,
	-0.0000000000000000000000017323564459105165Q,
	-0.0000000000000000000000000236061902449928Q,
	 0.0000000000000000000000000186498294171728Q,
	-0.0000000000000000000000000022180956242072Q,
	 0.0000000000000000000000000001297781974948Q,
	 0.0000000000000000000000000000011806974748Q,
	-0.0000000000000000000000000000011245843493Q,
	 0.0000000000000000000000000000001277085176Q,
	-0.0000000000000000000000000000000073914512Q,
	 0.0000000000000000000000000000000000113476Q,
	 0.0000000000000000000000000000000000463914Q,
	-0.0000000000000000000000000000000000053474Q,
	 0.0000000000000000000000000000000000003208Q,
	-0.0000000000000000000000000000000000000044Q,
	-0.0000000000000000000000000000000000000013Q,
	 0.0000000000000000000000000000000000000002Q,
      }};
      const auto s_eps = emsr::epsilon(std::real(a));
      auto ak = _Tp{1};
      auto gam = _Tp{0};
      for (auto k = 1u; k < s_c.size(); ++k)
	{
	  ak *= a;
	  auto term = s_c[k] * ak;
	  gam += term;
	  if (std::abs(term) < s_eps)
	    break;
	}
      return gam;
    }

  /**
   * Return the reciprocal of the Gamma function by infinite product.
   * The reciprocal of the Gamma function is given by
   * @f[
   *   \frac{1}{\Gamma(a)} = ae^{\gamma_E a}\Pi_{k=1}^{\infty}
   *                     (\left 1+\frac{a}{k}\right)e^{-a/k}
   * @f]
   */
  template<typename _Tp>
    _Tp
    gamma_reciprocal_prod(_Tp a)
    {
      const auto s_eps = emsr::epsilon(std::real(a));
      const auto s_gamma_e = emsr::egamma_v<_Tp>;
      const auto s_max_iter = 10000;
      auto gam = a * std::exp(s_gamma_e * a);
      for (auto k = 1u; k < s_max_iter; ++k)
	{
	  const auto rat = a / _Tp(k);
	  gam *= (_Tp{1} + rat) * std::exp(-rat);
	  if (std::abs(rat) < s_eps)
	    break;
	}
      return gam;
    }

  /**
   * Return the reciprocal of the Gamma function:
   * @f[
   *   \frac{1}{\Gamma(a)}
   * @f]
   *
   * @param a The argument of the reciprocal of the gamma function.
   * @return  The reciprocal of the gamma function.
   */
  template<typename _Tp>
    _Tp
    gamma_reciprocal(_Tp a)
    {
      using _Real = emsr::num_traits_t<_Tp>;

      if (std::isnan(a))
	return emsr::quiet_NaN(a);
      else
	{
	  const auto s_pi = emsr::pi_v<_Tp>;
	  const auto an = emsr::fp_is_integer(a);
	  if (an)
	    {
	      auto n = an();
	      if (n <= 0)
		return _Tp{0};
	      else if (n < int(emsr::detail::s_num_factorials<_Real>))
		return _Tp{1}
		    / _Real(emsr::detail::s_factorial_table[n - 1].factorial);
	      else
	        {
		  auto k = int(emsr::detail::s_num_factorials<_Real>);
		  auto rgam = _Tp{1}
			      / _Real(emsr::detail::s_factorial_table[k - 1].factorial);
		  while (k < n && std::abs(rgam) > _Real{0})
		    rgam /= _Real(k++);
		  return rgam;
		}
	    }
	  else if (std::real(a) > _Real{1})
	    {
	      auto n = int(std::floor(std::real(a)));
	      auto nu = a - _Tp(n);
	      auto rgam = gamma_reciprocal_series(nu);
	      while (std::real(a) > _Real{1} && std::abs(rgam) > _Tp{0})
	        rgam /= (a -= _Real{1});
	      return rgam;
	    }
	  else if (std::real(a) > _Real{0})
	    return gamma_reciprocal_series(a);
	  else
	    return emsr::detail::sin_pi(a)
		 * emsr::detail::gamma(_Tp{1} - a) / s_pi;
	}
    }

  /**
   * @brief A structure for the gamma functions required by the Temme series
   * 	    expansions of @f$ N_\nu(x) @f$ and @f$ K_\nu(x) @f$.
   * @f[
   *   \Gamma_1 = \frac{1}{2\mu}
   * 	 \left[\frac{1}{\Gamma(1 - \mu)} - \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * and
   * @f[
   *   \Gamma_2 = \frac{1}{2}
   *     \left[\frac{1}{\Gamma(1 - \mu)} + \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * where @f$ -1/2 <= \mu <= 1/2 @f$ is @f$ \mu = \nu - N @f$ and @f$ N @f$.
   * is the nearest integer to @f$ \nu @f$.
   * The values of @f$ \Gamma(1 + \mu) @f$ and @f$ \Gamma(1 - \mu) @f$
   * are returned as well.
   *
   * The accuracy requirements on this are high for @f$ |\mu| < 0 @f$.
   */
  template<typename _Tp>
    struct gamma_temme_t
    {
      /// The input parameter of the gamma functions
      _Tp mu_arg;
      /// The output function @f$ 1/\Gamma(1 + \mu) @f$
      _Tp gamma_plus_value;
      /// The output function @f$ 1/\Gamma(1 - \mu) @f$
      _Tp gamma_minus_value;
      /// The output function @f$ \Gamma_1(\mu) @f$
      _Tp gamma_1_value;
      /// The output function @f$ \Gamma_2(\mu) @f$
      _Tp gamma_2_value;
    };

  /**
   * @brief Compute the gamma functions required by the Temme series
   * 	    expansions of @f$ N_\nu(x) @f$ and @f$ K_\nu(x) @f$.
   * @f[
   *   \Gamma_1 = \frac{1}{2\mu}
   * 	 \left[\frac{1}{\Gamma(1 - \mu)} - \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * and
   * @f[
   *   \Gamma_2 = \frac{1}{2}
   *     \left[\frac{1}{\Gamma(1 - \mu)} + \frac{1}{\Gamma(1 + \mu)}\right]
   * @f]
   * where @f$ -1/2 <= \mu <= 1/2 @f$ is @f$ \mu = \nu - N @f$ and @f$ N @f$.
   * is the nearest integer to @f$ \nu @f$.
   * The values of @f$ \Gamma(1 + \mu) @f$ and @f$ \Gamma(1 - \mu) @f$
   * are returned as well.
   *
   * The accuracy requirements on this are exquisite.
   *
   * @param mu     The input parameter of the gamma functions.
   * @param[out] gam1   The output function @f$ \Gamma_1(\mu) @f$
   * @param[out] gam2   The output function @f$ \Gamma_2(\mu) @f$
   * @param[out] gamp  The output function @f$ 1/\Gamma(1 + \mu) @f$
   * @param[out] gamm  The output function @f$ 1/\Gamma(1 - \mu) @f$
   */
  template<typename _Tp>
    emsr::gamma_temme_t<_Tp>
    gamma_temme(_Tp mu)
    {
      using gammat_t = emsr::gamma_temme_t<_Tp>;
      const auto s_eps = emsr::epsilon(mu);
      const auto s_gamma_E = emsr::egamma_v<_Tp>;

      if (std::abs(mu) < s_eps)
	return gammat_t{mu, _Tp{1}, _Tp{1}, -s_gamma_E, _Tp{1}};
      else
	{
	  _Tp gamp, gamm;
	  if (std::real(mu) <= _Tp{0})
	    {
	      gamp = gamma_reciprocal_series(_Tp{1} + mu);
	      gamm = -gamma_reciprocal_series(-mu) / mu;
	    }
	  else
	    {
	      gamp = gamma_reciprocal_series(mu) / mu;
	      gamm = gamma_reciprocal_series(_Tp{1} - mu);
	    }
	  auto gam1 = (gamm - gamp) / (_Tp{2} * mu);
	  auto gam2 = (gamm + gamp) / _Tp{2};
	  return gammat_t{mu, gamp, gamm, gam1, gam2};
	}
    }

  template<typename _Tp>
    gamma_temme_t<_Tp>
    gamma_temme_std(_Tp mu)
    {
      const auto s_eps = emsr::epsilon(mu);
      const auto s_gamma_E = emsr::egamma_v<_Tp>;
      auto gamp = _Tp{1} / std::tgamma(_Tp{1} + mu);
      auto gamm = _Tp{1} / std::tgamma(_Tp{1} - mu);
      auto gam1 = (std::abs(mu) < s_eps)
		  ? -s_gamma_E
		  : (gamm - gamp) / (_Tp{2} * mu);
      auto gam2 = (gamm + gamp) / _Tp{2};

      return gamma_temme_t<_Tp>{mu, gamp, gamm, gam1, gam2};
    }


template<typename _Tp>
  void
  plot_gamma_reciprocal(_Tp proto)
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();
    const auto del = _Tp{1} / _Tp{100};
    for (auto k = -500; k <= 1000; ++k)
      {
	auto a = k * del;
	auto gammar = gamma_reciprocal(a);
	std::cout << ' ' << std::setw(width) << a
		<< ' ' << std::setw(width) << gammar
    		<< '\n';
      }
    std::cout << "\n\n";
  }


template<typename _Tp>
  void
  test_gamma_reciprocal(_Tp proto)
  {
    using _Val = _Tp;
    using _Real = emsr::num_traits_t<_Val>;

    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::size_t n = 50;
    auto c = gamma_reciprocal_series_coef<_Real>(n, proto);

    std::cout << '\n'
	      << ' ' << std::setw(4) << "k"
	      << ' ' << std::setw(width) << "c"
    	      << '\n';
    for (auto k = 0u; k < c.size(); ++k)
      std::cout << ' ' << std::setw(4) << k
		<< ' ' << std::setw(width) << c[k]
    		<< '\n';

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "1/G(a) ser"
	      << ' ' << std::setw(width) << "1/G(a) prd"
	      << ' ' << std::setw(width) << "1/std::tgm"
	      << ' ' << std::setw(width) << "1/G(a)"
	      << ' ' << std::setw(width) << "del ser"
	      << ' ' << std::setw(width) << "del prd"
	      << ' ' << std::setw(width) << "del"
    	      << '\n';
    const auto del = _Tp{1} / _Tp{100};
    for (auto k = -500; k <= 1000; ++k)
      {
	auto a = k * del;
	auto gammargs = gamma_reciprocal_series(a);
	auto gammargp = gamma_reciprocal_prod(a);
	auto gammarstd = _Tp{1} / std::tgamma(a);
	auto gammar = gamma_reciprocal(a);
	std::cout << ' ' << std::setw(width) << a
		<< ' ' << std::setw(width) << gammargs
		<< ' ' << std::setw(width) << gammargp
		<< ' ' << std::setw(width) << gammarstd
		<< ' ' << std::setw(width) << gammar
		<< ' ' << std::setw(width) << (gammargs - gammarstd) / gammarstd
		<< ' ' << std::setw(width) << (gammargp - gammarstd) / gammarstd
		<< ' ' << std::setw(width) << (gammar - gammarstd) / gammarstd
    		<< '\n';
      }
  }

template<typename _Tp>
  void
  test_gamma_temme(_Tp proto)
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "1/Gamma(1+mu)"
	      << ' ' << std::setw(width) << "1/Gamma(1-mu)"
	      << ' ' << std::setw(width) << "Gamma_1"
	      << ' ' << std::setw(width) << "Gamma_2"
	      << ' ' << std::setw(width) << "1/std::tgamma(1+mu)"
	      << ' ' << std::setw(width) << "1/std::tgamma(1-mu)"
	      << ' ' << std::setw(width) << "tgamma_1"
	      << ' ' << std::setw(width) << "tgamma_2"
	      << ' ' << std::setw(width) << "delta (1+mu)"
	      << ' ' << std::setw(width) << "delta (1-mu)"
	      << ' ' << std::setw(width) << "delta 1"
	      << ' ' << std::setw(width) << "delta 2"
    	      << '\n';
    const auto del = _Tp{1} / _Tp{100};
    for (auto k = -100; k <= 100; ++k)
      {
	auto mu = k * del;
	auto tggnu = gamma_temme(mu);
	auto tgstd = gamma_temme_std(mu);
	std::cout << ' ' << std::setw(width) << mu
		<< ' ' << std::setw(width) << tggnu.gamma_plus_value
		<< ' ' << std::setw(width) << tggnu.gamma_minus_value
		<< ' ' << std::setw(width) << tggnu.gamma_1_value
		<< ' ' << std::setw(width) << tggnu.gamma_2_value
		<< ' ' << std::setw(width) << tgstd.gamma_plus_value
		<< ' ' << std::setw(width) << tgstd.gamma_minus_value
		<< ' ' << std::setw(width) << tgstd.gamma_1_value
		<< ' ' << std::setw(width) << tgstd.gamma_2_value
		<< ' ' << std::setw(width) << (tggnu.gamma_plus_value - tgstd.gamma_plus_value) / tgstd.gamma_plus_value
		<< ' ' << std::setw(width) << (tggnu.gamma_minus_value - tgstd.gamma_minus_value) / tgstd.gamma_minus_value
		<< ' ' << std::setw(width) << (tggnu.gamma_1_value - tgstd.gamma_1_value) / tgstd.gamma_1_value
		<< ' ' << std::setw(width) << (tggnu.gamma_2_value - tgstd.gamma_2_value) / tgstd.gamma_2_value
    		<< '\n';
      }
  }

int
main()
{
  plot_gamma_reciprocal(1.0);

  test_gamma_reciprocal(1.0f);

  test_gamma_reciprocal(1.0);

  test_gamma_reciprocal(1.0l);

  //test_gamma_reciprocal(1.0q);

  //test_gamma_reciprocal(mpfr::mpreal(1,128));

  test_gamma_temme(1.0f);

  test_gamma_temme(1.0);

  test_gamma_temme(1.0l);

  //test_gamma_temme(1.0q);
}
