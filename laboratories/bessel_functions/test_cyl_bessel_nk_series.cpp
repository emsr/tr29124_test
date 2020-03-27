/**
 *
 */

#include <stdexcept>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <ext/math_constants.h>
#include <bits/numeric_limits.h>

  template<typename _Tp>
    struct __bessel_nk_series_t
    {
      _Tp _Z_mu;
      _Tp _Z_mup1;
    };

  template<typename _Tp>
    __gnu_cxx::__gamma_temme_t<_Tp>
    __gamma_temme(_Tp __mu)
    {
      using __gammat_t = __gnu_cxx::__gamma_temme_t<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__mu);
      const auto _S_gamma_E = __gnu_cxx::numbers::__gamma_e_v<_Tp>;

      if (std::abs(__mu) < _S_eps)
	return __gammat_t{__mu, _Tp{1}, _Tp{1}, -_S_gamma_E, _Tp{1}};
      else
	{
	  _Tp __gamp, __gamm;
	  if (std::real(__mu) <= _Tp{0})
	    {
	      __gamp = std::__detail::__gamma_reciprocal_series(_Tp{1} + __mu);
	      __gamm = -std::__detail::__gamma_reciprocal_series(-__mu) / __mu;
	    }
	  else
	    {
	      __gamp = std::__detail::__gamma_reciprocal_series(__mu) / __mu;
	      __gamm = std::__detail::__gamma_reciprocal_series(_Tp{1} - __mu);
	    }
	  const auto __gam1 = (__gamm - __gamp) / (_Tp{2} * __mu);
	  const auto __gam2 = (__gamm + __gamp) / _Tp{2};
	  return __gammat_t{__mu, __gamp, __gamm, __gam1, __gam2};
	}
    }

  template<typename _Tp>
    __bessel_nk_series_t<_Tp>
    old_n(_Tp __nu, _Tp __x, int __max_iter = 10000)
    {
      const auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Tp>;
      const int __n = std::nearbyint(__nu);
      const auto __mu = __nu - _Tp(__n);
      const auto __mu2 = __mu * __mu;
      const auto __xi = _Tp{1} / __x;
      const auto __xi2 = _Tp{2} * __xi;
      const auto __x2 = __x / _Tp{2};
      const auto __pimu = _S_pi * __mu;
      const auto __fact = (std::abs(__pimu) < _S_eps
			? _Tp{1}
			: __pimu / std::sin(__pimu));
      auto __d = -std::log(__x2);
      auto __e = __mu * __d;
      const auto __fact2 = (std::abs(__e) < _S_eps
			 ? _Tp{1}
			 : std::sinh(__e) / __e);
      const auto __gamt = __gamma_temme(__mu);
      auto __ff = (_Tp{2} / _S_pi) * __fact
		* (__gamt.__gamma_1_value * std::cosh(__e)
		 + __gamt.__gamma_2_value * __fact2 * __d);
      __e = std::exp(__e);
      auto __p = __e / (_S_pi * __gamt.__gamma_plus_value);
      auto __q = _Tp{1} / (__e * _S_pi * __gamt.__gamma_minus_value);
      const auto __pimu2 = __pimu / _Tp{2};
      const auto __fact3 = (std::abs(__pimu2) < _S_eps
			 ? _Tp{1} : std::sin(__pimu2) / __pimu2 );
      const auto __r = _S_pi * __pimu2 * __fact3 * __fact3;
      auto __c = _Tp{1};
      __d = -__x2 * __x2;
      auto __sum = __ff + __r * __q;
      auto __sum1 = __p;
      int __i;
      for (__i = 1; __i <= __max_iter; ++__i)
	{
	  __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
	  __c *= __d / _Tp(__i);
	  __p /= _Tp(__i) - __mu;
	  __q /= _Tp(__i) + __mu;
	  const auto __del = __c * (__ff + __r * __q);
	  __sum += __del;
	  const auto __del1 = __c * __p - _Tp(__i) * __del;
	  __sum1 += __del1;
	  if (std::abs(__del) < _S_eps * (_Tp{1} + std::abs(__sum)))
	    break;
	}
      if (__i > __max_iter)
	std::__throw_runtime_error(__N("__cyl_bessel_nk_series: "
				       "Series failed to converge"));

      auto _Nmu = -__sum;
      auto _Nnu1 = -__sum1 * __xi2;

      return {_Nmu, _Nnu1};
    }

  template<typename _Tp>
    __bessel_nk_series_t<_Tp>
    old_k(_Tp __nu, _Tp __x, int __max_iter = 10000)
    {
      const auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Tp>;
      const int __n = std::nearbyint(__nu);
      const auto __mu = __nu - _Tp(__n);
      const auto __mu2 = __mu * __mu;
      const auto __xi = _Tp{1} / __x;
      const auto __xi2 = _Tp{2} * __xi;
      const auto __x2 = __x / _Tp{2};
      const auto __pimu = _S_pi * __mu;
      const auto __fact = (std::abs(__pimu) < _S_eps
			? _Tp{1}
			: __pimu / std::sin(__pimu));
      auto __d = -std::log(__x2);
      auto __e = __mu * __d;
      const auto __fact2 = (std::abs(__e) < _S_eps
			 ? _Tp{1}
			 : std::sinh(__e) / __e);
      const auto __gamt = __gamma_temme(__mu);
      auto __ff = __fact
		* (__gamt.__gamma_1_value * std::cosh(__e)
		 + __gamt.__gamma_2_value * __fact2 * __d);
      auto __sum = __ff;
      __e = std::exp(__e);
      auto __p = __e / (_Tp{2} * __gamt.__gamma_plus_value);
      auto __q = _Tp{1} / (_Tp{2} * __e * __gamt.__gamma_minus_value);
      auto __c = _Tp{1};
      __d = __x2 * __x2;
      auto __sum1 = __p;
      int __i;
      for (__i = 1; __i <= __max_iter; ++__i)
	{
	  __ff = (__i * __ff + __p + __q) / (__i * __i - __mu2);
	  __c *= __d / _Tp(__i);
	  __p /= _Tp(__i) - __mu;
	  __q /= _Tp(__i) + __mu;
	  const auto __del = __c * __ff;
	  __sum += __del;
	  const auto __del1 = __c * (__p - _Tp(__i) * __ff);
	  __sum1 += __del1;
	  if (std::abs(__del) < _S_eps * std::abs(__sum))
	    break;
	}
      if (__i > __max_iter)
	std::__throw_runtime_error(__N("__cyl_bessel_ik_steed: "
				       "K-series failed to converge"));
      auto _Kmu = __sum;
      auto _Kmu1 = __sum1 * __xi2;

      return {_Kmu, _Kmu1};
    }

  /**
   * This routine computes the dominant cylindrical bessel function solutions
   * by series summation for order @f$ |\mu| < 1/2 @f$.
   *
   * @param __mu The order of the Bessel functions @f$ |\mu| < 1/2 @f$.
   * @param __x  The argument of the Bessel functions.
   * @param __modified If true solve for the modified Bessel function
   *                   @f$ K_\mu @f$, otherwise solve for the Neumann function
   *                   @f$ N_\mu @f$,
   * @return A structure containing Z_{\mu} and Z_{\mu+1}.
   */
  template<typename _Tp>
    __bessel_nk_series_t<_Tp>
    __cyl_bessel_nk_series(_Tp __mu, _Tp __x, bool __modified = false,
			   int __max_iter = 100)
    {
      const auto _S_eps = __gnu_cxx::__epsilon<_Tp>();
      const auto _S_pi = __gnu_cxx::numbers::__pi_v<_Tp>;
      const auto __xi = _Tp{1} / __x;
      const auto __x2 = __x / _Tp{2};

      const auto __fact = _Tp{1} / std::__detail::__sinc_pi(__mu);
      const auto __lx2 = -std::log(__x2);
      const auto __arg = __mu * __lx2;
      const auto __fact2 = std::__detail::__sinhc(__arg);
      const auto __gamt = std::__detail::__gamma_temme(__mu);
      const auto __norm = __modified ? _Tp{-1} : _Tp{2} / _S_pi;
      auto __ff = __norm * __fact
		* (__gamt.__gamma_1_value * std::cosh(__arg)
		 + __gamt.__gamma_2_value * __fact2 * __lx2);
      const auto __e = std::exp(__arg);
      auto __p = __norm * __e / (_Tp{2} * __gamt.__gamma_plus_value);
      auto __q = __norm / (__e * _Tp{2} * __gamt.__gamma_minus_value);
      const auto __fact3 = __modified
			 ? _Tp{0}
			 : std::__detail::__sinc_pi(__mu / _Tp{2});
      const auto __r = __modified
		     ? _Tp{0}
		     : __fact3 * __fact3 * _S_pi * _S_pi * __mu / _Tp{2};
      auto __c = _Tp{1};
      const auto __d = __modified ? __x2 * __x2 : -__x2 * __x2;
      auto __sum_mu = __ff + __r * __q;
      auto __sum_mup1 = __p;
      int __i;
      for (__i = 1; __i <= __max_iter; ++__i)
	{
	  __ff = (__i * __ff + __p + __q)
	       / ((_Tp(__i) - __mu) * (_Tp(__i) + __mu));
	  __c *= __d / _Tp(__i);
	  __p /= _Tp(__i) - __mu;
	  __q /= _Tp(__i) + __mu;
	  const auto __del_mu = __c * (__ff + __r * __q);
	  __sum_mu += __del_mu;
	  const auto __del_mup1 = __c * __p - _Tp(__i) * __del_mu;
	  __sum_mup1 += __del_mup1;
	  if (std::abs(__del_mu) < _S_eps * std::abs(__sum_mu))
	    break;
	}
      if (__i > __max_iter)
	std::__throw_runtime_error(__N("__cyl_bessel_nk_series: "
				       "Series failed to converge"));
      auto _N_mu = -__sum_mu;
      auto _N_mup1 = -_Tp{2} * __xi * __sum_mup1;

      return {_N_mu, _N_mup1};
    }

template<typename _Tp>
  void
  test_cyl_bessel_nk_series()
  {
    const auto p = std::numeric_limits<_Tp>::digits10;
    std::cout.precision(p);
    const auto w = 8 + std::cout.precision();

    std::cout << "\ncyl_neumann\n";
    std::cout << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "new N_mu"
	      << ' ' << std::setw(w) << "new N_mup1"
	      << ' ' << std::setw(w) << "old N_mu"
	      << ' ' << std::setw(w) << "old N_mup1"
	      << ' ' << std::setw(w) << "new N_mup1 - std"
	      << ' ' << std::setw(w) << "N_mu (new-lab)/lab"
	      << ' ' << std::setw(w) << "N'_mu (new-lab)/lab"
	      << '\n';
    for (auto nu : {_Tp{0}, _Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}})
      {
	for (int i = 1; i < 20; ++i)
	  {
	    const auto x = i * 0.1;
	    const auto jn = std::__detail::__cyl_bessel_jn(nu, x);
	    const auto np1 = std::cyl_neumann(nu + 1, x);
	    const auto nn = __cyl_bessel_nk_series(nu, x);
	    const auto Np_mu = nu * nn._Z_mu / x - nn._Z_mup1;
	    const auto no = old_n(nu, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << nn._Z_mu
		      << ' ' << std::setw(w) << nn._Z_mup1
		      << ' ' << std::setw(w) << no._Z_mu
		      << ' ' << std::setw(w) << no._Z_mup1
		      << ' ' << std::setw(w) << nn._Z_mup1 - np1
		      << ' ' << std::setw(w) << (nn._Z_mu - jn.__N_value) / jn.__N_value
		      << ' ' << std::setw(w) << (Np_mu - jn.__N_deriv) / jn.__N_deriv
		      << '\n';
	  }
      }

    std::cout << "\ncyl_bessel_k\n";
    std::cout << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "new K_mu"
	      << ' ' << std::setw(w) << "new K_mup1"
	      << ' ' << std::setw(w) << "old K_mu"
	      << ' ' << std::setw(w) << "old K_mup1"
	      << ' ' << std::setw(w) << "new K_mup1 - std"
	      << ' ' << std::setw(w) << "K_mu (new-lab)/lab"
	      << ' ' << std::setw(w) << "K'_mu (new-lab)/lab"
	      << '\n';
    for (auto nu : {_Tp{0}, _Tp{1}/_Tp{3}, _Tp{1}/_Tp{2}})
      {
	for (int i = 1; i < 20; ++i)
	  {
	    const auto x = i * 0.1;
	    const auto ik = std::__detail::__cyl_bessel_ik(nu, x);
	    const auto kp1 = std::cyl_bessel_k(nu + 1, x);
	    const auto kn = __cyl_bessel_nk_series(nu, x, true);
	    const auto Np_mu = nu * kn._Z_mu / x - kn._Z_mup1;
	    const auto ko = old_k(nu, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << kn._Z_mu
		      << ' ' << std::setw(w) << kn._Z_mup1
		      << ' ' << std::setw(w) << ko._Z_mu
		      << ' ' << std::setw(w) << ko._Z_mup1
		      << ' ' << std::setw(w) << kn._Z_mup1 - kp1
		      << ' ' << std::setw(w) << (kn._Z_mu - ik.__K_value) / ik.__K_value
		      << ' ' << std::setw(w) << (Np_mu - ik.__K_deriv) / ik.__K_deriv
		      << '\n';
	  }
      }
  }


int
main()
{
  test_cyl_bessel_nk_series<double>();

  return 0;
}
