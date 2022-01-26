/**
 *
 */

#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

#include <statistics.h>

#include <emsr/summation.h>
#include <emsr/fp_type_util.h>
#include <emsr/sf_zeta.h>

#include <3rdparty/lerchphi/Source/lerchphi.h>

  /**
   * A functor for a vanWijnGaarden compressor.
   * vanWijnGaarden requires:
   *   _Tp operator()(int) that returns a term in the original defining series.
   */
  template<typename _Tp>
    class lerch_term
    {
    public:

      using value_type = _Tp;

      lerch_term(value_type z, value_type s, value_type a)
      : _M_z{z}, _M_s{s}, _M_a{a}
      { }

      value_type
      operator()(std::size_t i) const
      {
	return std::pow(_M_z, value_type(i))
	     / std::pow(_M_a + value_type(i), _M_s);
      }

    private:

      value_type _M_z;
      value_type _M_s;
      value_type _M_a;
    };

  /**
   * This function blows up on nonpositive integeral a.
   */
  template<typename _Tp>
    _Tp
    lerch_sum(_Tp z, _Tp s, _Tp a)
    {
      const auto s_nan = emsr::quiet_NaN(s);
      const auto s_eps = emsr::epsilon(s);

      const auto aint = emsr::fp_is_integer(a);
      if (aint && aint() <= 0)
	return s_nan;
      else if (std::abs(std::abs(z) - _Tp{1}) < s_eps
		&& std::real(s) <= _Tp{1} + s_eps)
	return s_nan;
      else if (std::abs(z) > _Tp{1} + s_eps)
	return s_nan;
      else
	{
	  constexpr auto s_maxit = 100000u;
	  auto zpow = _Tp{1};
	  auto sum = std::pow(a, -s);
	  for (auto k = 1u; k < s_maxit; ++k)
	    {
	      zpow *= z;
	      auto term = zpow * std::pow(a + k, -s);
	      sum += term;
	      if (std::abs(term / sum) < s_eps)
		break;
	    }
	  return sum;
	}
    }

  /**
   * This function blows up on nonpositive integeral a.
   */
  template<typename _Tp>
    _Tp
    lerch_vanwijngaarden_sum(_Tp z, _Tp s, _Tp a)
    {
      const auto s_nan = emsr::quiet_NaN(s);
      const auto s_eps = emsr::epsilon(s);

      const auto aint = emsr::fp_is_integer(a);
      if (aint && aint() <= 0)
	return s_nan;
      else if (std::abs(std::abs(z) - _Tp{1}) < s_eps
		&& std::real(s) <= _Tp{1} + s_eps)
	return s_nan;
      else if (std::abs(z) > _Tp{1} + s_eps)
	return s_nan;
      else if (z < _Tp{0})
	{
	  constexpr auto s_maxit = 100000u;
	  using lerch_t = lerch_term<_Tp>;
	  auto lerch_fun = lerch_t(z, s, a);
	  emsr::VanWijngaardenSum<_Tp> sum;
	  for (auto k = 0u; k < s_maxit; ++k)
	    {
	      auto temp = lerch_fun(k);
	      sum += temp;
	      if (std::abs(temp / sum) < s_eps)
		break;
	    }
	  return sum();
	}
      else
	{
	  constexpr auto s_maxit = 100000u;
	  using lerch_t = lerch_term<_Tp>;
	  auto lerch_fun = lerch_t(z, s, a);
	  emsr::VanWijngaardenCompressor<lerch_t> term(lerch_fun);
	  emsr::VanWijngaardenSum<_Tp> sum;
	  for (auto k = 0u; k < s_maxit; ++k)
	    {
	      auto temp = term[k];
	      sum += temp;
	      if (std::abs(temp / sum) < s_eps)
		break;
	    }
	  return sum();
	}
    }

  /**
   * This function blows up on nonpositive integeral a.
   *  As usual, the binomial coefficient kills this for practical purposes.
   */
  template<typename _Tp>
    _Tp
    lerch_double_sum(_Tp z, _Tp s, _Tp a)
    {
      const auto s_nan = emsr::quiet_NaN(s);
      const auto s_eps = emsr::epsilon(s);

      const auto aint = emsr::fp_is_integer(a);
      if (aint && aint() <= 0)
	return s_nan;
      else if (std::abs(std::abs(z) - _Tp{1}) < s_eps
		&& std::real(s) <= _Tp{1} + s_eps)
	return s_nan;
      else if (std::abs(z) > _Tp{1} + s_eps)
	return s_nan;
      else
	{
	  constexpr auto s_maxit = 10000u;
	  auto lerch = std::pow(a, -s);
	  const auto zfrac = -z / (_Tp{1} - z);
	  auto zfact = _Tp{1};
	  for (auto n = 1u; n < s_maxit; ++n)
	    {
	      auto term = std::pow(a, -s);
	      auto binomial = _Tp{1};
	      emsr::VanWijngaardenSum<_Tp> sum(term);
	      for (auto k = 1; k <= n; ++k)
		{
		  binomial *= -_Tp(n - k + 1) / _Tp(k);
		  term *= z * binomial * std::pow(a + k, -s);
		  sum += term;
		}
	      zfact *= zfrac;
	      lerch += zfact * sum();
	      if (std::abs(zfact * sum() / lerch) < s_eps)
		break;
	    }
	  lerch /= (_Tp{1} - z);
	  return lerch;
	}
    }

  /**
   * Try the WenigerDelta<MonotoneVanWijngaarden> composition.
   */
  template<typename _Tp>
    _Tp
    lerch_delta_vanwijngaarden_sum(_Tp z, _Tp s, _Tp a)
    {
      const auto s_eps = emsr::epsilon(s);
      constexpr auto s_maxit = 1000u;

      emsr::WenigerDeltaSum<emsr::VanWijngaardenSum<_Tp>> _WDvW;
      if (z >= _Tp{0})
	{
	  using lerch_t = lerch_term<_Tp>;
	  using lerch_cmp_t = emsr::VanWijngaardenCompressor<lerch_t>;
	  auto _VwT = lerch_cmp_t(lerch_t(z, s, a));
	  for (auto k = 0u; k < s_maxit; ++k)
	    {
	      auto term = _VwT[k];
	      _WDvW += term;
	      if (std::abs(term) < s_eps * std::abs(_WDvW()))
		break;
	    }
	  return _WDvW();
	}
      else
	{
	  auto _LT = lerch_term<_Tp>(z, s, a);
	  for (auto k = 0u; k < s_maxit; ++k)
	    {
	      auto term = _LT(k);
	      _WDvW += term;
	      if (std::abs(term) < s_eps * std::abs(_WDvW()))
		break;
	    }
	  return _WDvW();
	}
    }

  /**
   * This function blows up on nonpositive integeral a.
   */
  template<typename _Tp>
    _Tp
    lerch_phi(_Tp z, _Tp s, _Tp a)
    {
      const auto s_nan = emsr::quiet_NaN(s);
      const auto s_eps = emsr::epsilon(s);

      if (std::isnan(z) || std::isnan(s) || std::isnan(a))
	return s_nan;
      else if (std::abs(std::abs(z) - _Tp{1}) < s_eps
		&& std::real(s) <= _Tp{1} + s_eps)
	return s_nan;
      else if (std::abs(z) > _Tp{1} + s_eps)
	return s_nan;
      else
	{
	  const auto aint = emsr::fp_is_integer(a);

	  const auto sint = emsr::fp_is_integer(s);
	  const bool tinyz = std::abs(z) < s_eps; // s_min?
	  const bool smallz = !tinyz && (std::abs(z) < _Tp{0.5});

	  if (aint && aint() <= 0)
	    return s_nan;
	  else if (a < _Tp{0})
	    {
	      if (sint)
		{
		  int sign = sint() % 2 == 0 ? +1 : -1;
		  if (tinyz)
		    return sign * _Tp{1} / std::pow(std::abs(a), s);
		  else
		    {
		      auto m = -int(std::floor(a));
		      auto a1 = a + _Tp(m);
		      auto sum1 = _Tp{0};
		      for (int i = 0; i < m; ++i)
			{
			  sum1 += sign * std::pow(std::abs(z), i)
				  / std::pow(std::abs(a + i), _Tp(sint()));
			  if (z < _Tp{0})
			    sign = -sign;
			}
		      auto sum = _Tp{0};
		      if (smallz)
			sum = lerch_sum(z, s, a1);
		      else
			sum
			  = lerch_delta_vanwijngaarden_sum(z, s, a1);
		      sign = 1;
		      if (z < _Tp{0} && m % 2 != 0)
			sign = -1;
		      return sum1
			   + sum * sign * std::pow(std::abs(z), m);
		    }
		}
	      else // s is not an integer - Phi is complex.
		return s_nan;
	    }
	  else if (tinyz)
	    return _Tp{1} / std::pow(a, s);
	  else // a > 0
	    {
	      if (smallz)
		return lerch_sum(z, s, a);
	      else
		return lerch_delta_vanwijngaarden_sum(z, s, a);
	    }
	}
    }

  /**
   * Return the Hurwitz zeta function by evaluating the Lerch trancendent:
   * @f[
   *   \zeta(s,a) = \Phi(1,s,a)
   * @f]
   * @param[in] s The argument s > 1
   * @param[in] a The parameter
   */
  template<typename _Tp>
    _Tp
    hurwitz_zeta_lerch(_Tp s, _Tp a)
    {
      return lerch_phi(_Tp{1}, s, a);
    }

  /**
   * Return the Riemann zeta function by evaluating the Lerch trancendent:
   * @f[
   *   \zeta(s) = \Phi(1,s,1)
   * @f]
   * @param[in] s The argument s > 1
   */
  template<typename _Tp>
    _Tp
    riemann_zeta_lerch(_Tp s)
    {
      return lerch_phi(_Tp{1}, s, _Tp{1});
    }

  /**
   * Return the Dirichlet beta function by evaluating the Lerch trancendent:
   * @f[
   *   \beta(s) = \frac{1}{2^s}\Phi(-1,s,\frac{1}{2})
   * @f]
   * @param[in] s The argument s > 1
   */
  template<typename _Tp>
    _Tp
    dirichlet_beta_lerch(_Tp s)
    {
      return lerch_phi(_Tp{-1}, s, _Tp{0.5L}) / std::pow(_Tp{2}, s);
    }

  /**
   * Return the Dirichlet eta function by evaluating the Lerch trancendent:
   * @f[
   *   \eta(s) = \Phi(-1,s,1)
   * @f]
   * @param[in] s The argument s > 1
   */
  template<typename _Tp>
    _Tp
    dirichlet_eta_lerch(_Tp s)
    {
      return lerch_phi(_Tp{-1}, s, _Tp{1});
    }

  /**
   * Return the Dirichlet lambda function by evaluating the Lerch trancendent:
   * @f[
   *   \beta(s) = \frac{1}{2^s}\Phi(1,s,\frac{1}{2})
   * @f]
   * @param[in] s The argument s > 1
   */
  template<typename _Tp>
    _Tp
    dirichlet_lambda_lerch(_Tp s)
    {
      return lerch_phi(_Tp{1}, s, _Tp{0.5L}) / std::pow(_Tp{2}, s);
    }

  /**
   * Return the polylog function by evaluating the Lerch trancendent:
   * @f[
   *   \L_s(z) = \Phi(z,s,1)
   * @f]
   * @param[in] s The argument s > 1
   */
  template<typename _Tp>
    _Tp
    polylog_lerch(_Tp s, _Tp z)
    {
      return lerch_phi(z, s, _Tp{1});
    }

  /**
   * Return the Legendre chi function by evaluating the Lerch trancendent:
   * @f[
   *   \chi_\nu(z) = \frac{z}{2^\nu}\Phi(z^2,\nu,1)
   * @f]
   * @param[in] s The argument s > 1
   */
  template<typename _Tp>
    _Tp
    legendre_chi(_Tp nu, _Tp z)
    {
      return z * lerch_phi(z * z, nu, _Tp{0.5L})
	   / std::pow(_Tp{2}, nu);
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    fermi_dirac_lerch(_Tp s, _Tp mu)
    {
      auto expmu = std::exp(-mu);
      auto gamsp1 = std::tgamma(_Tp{1} + s);
      return gamsp1 * lerch_phi(-expmu, _Tp{1} + s, _Tp{1}) / expmu;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    bose_einstein_lerch(_Tp s, _Tp mu)
    {
      auto expmu = std::exp(mu);
      auto gamsp1 = std::tgamma(_Tp{1} + s);
      return expmu * gamsp1 * lerch_phi(expmu, _Tp{1} + s, _Tp{1});
    }

  float
  lerch_phif(float z, float s, float a)
  { return lerch_phi<float>(z, s, a); }

  long double
  lerch_phil(long double z, long double s, long double a)
  { return lerch_phi<long double>(z, s, a); }

  template<typename _Tpz, typename _Tps, typename _Tpa>
    emsr::fp_promote_t<_Tpz, _Tps, _Tpa>
    lerch_phi(_Tpz z, _Tps s, _Tpa a)
    {
      using type = emsr::fp_promote_t<_Tpz, _Tps, _Tpa>;
      return lerch_phi<type>(z, s, a);
    }


struct lerch_testcase
{
  double phi;
  double z;
  double s;
  double a;
  int flag;
};

lerch_testcase
lerch_tests[12]
{
  { 1.0000000000000000e+00, -1.0000000000000000e+00,  2.0000000000000000e+00,  1.0000000000000000e+00, 1},
  { 1.0000000000000000e+00,  9.9999000000000005e-01,  2.0000000000000000e+00, -1.0000000000000000e+00, 2},
  { 1.0000000000000000e+00,  9.9999000000000005e-01,  2.2999999999999998e+00, -1.5000000000000000e+00, 3},
  { 1.8420680923134405e+01,  9.9999998999999995e-01,  1.0000000000000000e+00,  1.0000000000000000e+00, 0},
  { 1.6448253852467796e+00,  9.9999000000000005e-01,  2.0000000000000000e+00,  1.0000000000000000e+00, 0},
  { 8.2246832662591696e-01, -9.9999000000000005e-01,  2.0000000000000000e+00,  1.0000000000000000e+00, 0},
  { 9.5971489709979654e-04,  9.9999000000000005e-01,  2.0000000000000000e+00,  1.0000000000000000e+03, 0},
  { 1.4275808137603091e-01,  2.9999999999999999e-01,  2.0000000000000000e+00, -4.5000000000000000e+00, 0},
  { 1.0000025000111110e+00,  1.0000000000000001e-05,  2.0000000000000000e+00,  1.0000000000000000e+00, 0},
  { 9.9998425044098438e-01, -6.3000000000000000e-05,  2.0000000000000000e+00,  1.0000000000000000e+00, 0},
  { 6.5909228798196373e-01,  3.4709929976435479e-06,  1.0000000000000000e+00,  1.5172413793103448e+00, 0},
  { 2.5880201290103731e+17,  2.9999999999999997e-04,  2.0000000000000000e+00, -3.0000000000000102e+00, 0},
};

int
main()
{
  using Tp = double;
  constexpr auto s_nan = std::numeric_limits<Tp>::quiet_NaN();

  std::cout.precision(std::numeric_limits<Tp>::digits10);
  auto width = 8 + std::numeric_limits<Tp>::digits10;

  std::cout << "case " << std::setw(2) << "i"
	    << std::setw(width) << "z"
	    << std::setw(width) << "s"
	    << std::setw(width) << "a"
	    << std::setw(width) << "phi0"
	    << std::setw(6) << "flag"
	    << std::setw(width) << "phi"
	    << std::setw(width) << "phi-phi0"
	    << std::setw(width) << "lphi"
	    << std::setw(width) << "lphi-phi0"
	    << std::setw(width) << "phi-lphi"
	    << '\n';
  std::cout << "---- " << std::setw(2) << "-"
	    << std::setw(width) << "---------"
	    << std::setw(width) << "---------"
	    << std::setw(width) << "---------"
	    << std::setw(width) << "---------"
	    << std::setw(6) << "----"
	    << std::setw(width) << "---------"
	    << std::setw(width) << "---------"
	    << std::setw(width) << "---------"
	    << std::setw(width) << "---------"
	    << std::setw(width) << "---------"
	    << '\n';
  for (int i = 0; i < 12; ++i)
    {
      const auto& tcase = lerch_tests[i];
      std::cout << "case " << std::setw(2) << i + 1
                << std::setw(width) << tcase.z
                << std::setw(width) << tcase.s
                << std::setw(width) << tcase.a
                << std::setw(width) << tcase.phi
		<< std::setw(6) << tcase.flag;
      auto phi = Tp{0};
      try
	{
	  phi = lerch_phi(tcase.z, tcase.s, tcase.a);
	  auto test0 = phi - tcase.phi;
	  std::cout << std::setw(width) << phi
                    << std::setw(width) << test0;
	}
      catch (...)
	{
	  std::cout << std::setw(width) << "fail";
	  std::cout << std::setw(width) << "fail";
	  phi = s_nan;
	}
      double acc = 2 * std::numeric_limits<Tp>::epsilon();
      double lphi = 0.0;
      int iter = 0;
      auto ok = lerchphi(&tcase.z, &tcase.s, &tcase.a, &acc, &lphi, &iter);
      if (ok == 0)
	std::cout << std::setw(width) << lphi
		  << std::setw(width) << lphi - tcase.phi
		  << std::setw(width) << phi - lphi;
      else
	std::cout << std::setw(width) << "fail"
		  << std::setw(width) << "fail"
		  << std::setw(width) << "fail";
      if (ok != tcase.flag)
	std::cout << std::setw(12) << "flag error";

      std::cout << '\n';
    }

  auto s = 1.0;
  auto a = 2.0;
  std::cout << "\n Phi(z,1,2) Tests\n";
  std::cout << " s = " << std::setw(width) << s << '\n';
  std::cout << " a = " << std::setw(width) << a << '\n';
  for (int iz = -99; iz <= +99; ++iz)
    {
      auto z = 0.01 * iz;
      auto lerch1 = lerch_sum(z, s, a);
      auto lerch2 = lerch_vanwijngaarden_sum(z, s, a);
      //auto lerch3 = lerch_double_sum(z, s, a);
      auto lerch4 = lerch_delta_vanwijngaarden_sum(z, s, a);
      double acc = 2 * std::numeric_limits<Tp>::epsilon();
      double lphi = 0.0;
      int iter = 0;
      auto ok = lerchphi(&z, &s, &a, &acc, &lphi, &iter);
      if (ok != 0)
	lphi = s_nan;
      std::cout << ' ' << std::setw(width) << z
		<< ' ' << std::setw(width) << lerch1
		<< ' ' << std::setw(width) << lerch2
		//<< ' ' << std::setw(width) << lerch3
		<< ' ' << std::setw(width) << lerch4
		<< ' ' << std::setw(width) << lerch2 - lerch1
		<< ' ' << std::setw(width) << lerch4 - lerch1
		<< ' ' << std::setw(width) << lerch4 - lphi
		<< '\n';
    }

  s = 2.0;
  a = 1.0;
  std::cout << "\n Phi(z,2,1) Tests\n";
  std::cout << " s = " << std::setw(width) << s << '\n';
  std::cout << " a = " << std::setw(width) << a << '\n';
  for (int iz = -99; iz <= +99; ++iz)
    {
      auto z = 0.01 * iz;
      auto lerch1 = lerch_sum(z, s, a);
      auto lerch2 = lerch_vanwijngaarden_sum(z, s, a);
      //auto lerch3 = lerch_double_sum(z, s, a);
      auto lerch4 = lerch_delta_vanwijngaarden_sum(z, s, a);
      double acc = 2 * std::numeric_limits<Tp>::epsilon();
      double lphi = 0.0;
      int iter = 0;
      auto ok = lerchphi(&z, &s, &a, &acc, &lphi, &iter);
      if (ok != 0)
	lphi = s_nan;
      std::cout << ' ' << std::setw(width) << z
		<< ' ' << std::setw(width) << lerch1
		<< ' ' << std::setw(width) << lerch2
		//<< ' ' << std::setw(width) << lerch3
		<< ' ' << std::setw(width) << lerch4
		<< ' ' << std::setw(width) << lerch2 - lerch1
		<< ' ' << std::setw(width) << lerch4 - lerch1
		<< ' ' << std::setw(width) << lerch4 - lphi
		<< '\n';
    }

  std::cout << "\nDilogarithm Tests\n";
  _Statistics<Tp> dilog_stats;
  s = 2.0;
  a = 1.0;
  std::cout << '\n';
  std::cout << " s = " << std::setw(width) << s << '\n';
  std::cout << " a = " << std::setw(width) << a << '\n';
  std::cout << std::setw(width) << "z"
	    << std::setw(width) << "z Phi(z, s, 1)"
	    << std::setw(width) << "Li_s(z)"
	    << std::setw(width) << "zPhi - Li"
	    << '\n';
  for (int iz = -99; iz <= +99; ++iz)
    {
      auto z = 0.01 * iz;
      auto zlerch = s_nan;
      try
	{
	  zlerch = z * lerch_phi(z, s, a);
	}
      catch (...)
	{
	  std::cout << " fail\n";
	}
      auto dilog = emsr::dilog(z);
      auto delta = zlerch - dilog;
      dilog_stats << delta;
      std::cout << ' ' << std::setw(width) << z
		<< ' ' << std::setw(width) << zlerch
		<< ' ' << std::setw(width) << dilog
		<< ' ' << std::setw(width) << delta
		<< '\n';
    }
  std::cout << "// mean(Phi - zeta)    : " << dilog_stats.mean() << '\n';
  std::cout << "// variance(Phi - zeta): " << dilog_stats.variance() << '\n';
  std::cout << "// stddev(Phi - zeta)  : " << dilog_stats.std_deviation() << '\n';

  std::cout << "\nRiemann Zeta Tests\n";
  _Statistics<Tp> riemann_stats;
  auto z = 1.0;
  a = 1.0;
  std::cout << '\n';
  std::cout << " z = " << std::setw(width) << z << '\n';
  std::cout << " a = " << std::setw(width) << a << '\n';
  std::cout << std::setw(width) << "s"
	    << std::setw(width) << "Phi(1, s, 1)"
	    << std::setw(width) << "zeta(s)"
	    << std::setw(width) << "Phi - zeta"
	    << '\n';
  for (int is = 101; is <= 200; ++is)
    {
      auto s = 0.01 * is;
      auto lerch = s_nan;
      try
	{
	  lerch = lerch_phi(z, s, a);
	}
      catch (...)
	{
	}
      auto zeta = emsr::riemann_zeta(s);
      auto delta = lerch - zeta;
      riemann_stats << delta;
      std::cout << ' ' << std::setw(width) << s
		<< ' ' << std::setw(width) << lerch
		<< ' ' << std::setw(width) << zeta
		<< ' ' << std::setw(width) << delta
		<< '\n';
    }
  std::cout << "// mean(Phi - zeta)    : " << riemann_stats.mean() << '\n';
  std::cout << "// variance(Phi - zeta): " << riemann_stats.variance() << '\n';
  std::cout << "// stddev(Phi - zeta)  : " << riemann_stats.std_deviation() << '\n';

  std::cout << "\nHurwitz Zeta Tests\n";
  for (int ia = 1; ia <= 10; ++ia)
    {
      _Statistics<Tp> hurwitz_stats;
      auto a = 1.0 * ia;
      auto z = 1.0;
      std::cout << '\n';
      std::cout << " z = " << std::setw(width) << z << '\n';
      std::cout << " a = " << std::setw(width) << a << '\n';
      std::cout << std::setw(width) << "s"
		<< std::setw(width) << "Phi(1, s, a)"
		<< std::setw(width) << "zeta(s, a)"
		<< std::setw(width) << "Phi - zeta"
		<< '\n';
      for (int is = 101; is <= 200; ++is)
	{
	  auto s = 0.01 * is;
	  auto lerch = s_nan;
	  try
	    {
	      lerch = lerch_phi(z, s, a);
	    }
	  catch (...)
	    {
	    }
	  auto zeta = emsr::hurwitz_zeta(s, a);
	  auto delta = lerch - zeta;
	  hurwitz_stats << delta;
	  std::cout << ' ' << std::setw(width) << s
		    << ' ' << std::setw(width) << lerch
		    << ' ' << std::setw(width) << zeta
		    << ' ' << std::setw(width) << delta
		    << '\n';
	}
      std::cout << "// mean(Phi - zeta)    : " << hurwitz_stats.mean() << '\n';
      std::cout << "// variance(Phi - zeta): " << hurwitz_stats.variance() << '\n';
      std::cout << "// stddev(Phi - zeta)  : " << hurwitz_stats.std_deviation() << '\n';
    }

  std::cout << '\n';
  for (int ia = 1; ia <= 10; ++ia)
    {
      auto a = 1.0 * ia;
      std::cout << "\n a = " << std::setw(width) << a << '\n';
      for (int is = 0; is <= 50; ++is)
	{
	  auto s = 0.1 * is;
	  std::cout << "\n s = " << std::setw(width) << s << '\n' << '\n';
	  for (int iz = -99; iz <= +99; ++iz)
	    {
	      auto z = 0.01 * iz;
	      auto lerch1 = lerch_sum(z, s, a);
	      auto lerch2 = lerch_vanwijngaarden_sum(z, s, a);
	      //auto lerch3 = lerch_double_sum(z, s, a);
	      auto lerch4 = lerch_delta_vanwijngaarden_sum(z, s, a);
	      double acc = 2 * std::numeric_limits<Tp>::epsilon();
	      double lphi = 0.0;
	      int iter = 0;
	      auto ok = lerchphi(&z, &s, &a, &acc, &lphi, &iter);
	      if (ok != 0)
	        lphi = s_nan;
	      std::cout << ' ' << std::setw(width) << z
			<< ' ' << std::setw(width) << lerch1
			<< ' ' << std::setw(width) << lerch2
			//<< ' ' << std::setw(width) << lerch3
			<< ' ' << std::setw(width) << lerch4
			<< ' ' << std::setw(width) << lerch2 - lerch1
			<< ' ' << std::setw(width) << lerch4 - lerch1
			<< ' ' << std::setw(width) << lerch4 - lphi
			<< '\n';
	    }
	}
    }

  //auto lerch1 = lerch_vanwijngaarden_sum(-0.75, Tp{1}, Tp{2});
  //auto lerch2 = lerch_vanwijngaarden_sum(-0.5, Tp{0}, Tp{1});
}
