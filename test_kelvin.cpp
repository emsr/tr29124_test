/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_kelvin test_kelvin.cpp
./test_kelvin > test_kelvin.txt
*/

#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <bits/summation.h>
#include <bits/sf_gamma.tcc>
#include <ext/math_util.h>
#include <bits/complex_util.h>

namespace std
{
namespace __detail
{

  /**
   * Data structure for Kelvin functions.
   */
  template<typename _Real>
    struct _KelvinState
    {
      _Real __nu;
      _Real __x;
      _Real __ber;
      _Real __bei;
      _Real __ker;
      _Real __kei;
    };


  /**
   * Compute the Kelvin functions by series summation:
   * @f[
   *    ber(x) = \sum_{k=0}^\infty \frac{-(x/4)^4}{[(1/2)_k]^2[(1)_k]^2}
   * @f]
   * @f[
   *    bei(x) = \frac{x^2}{4}
   *             \sum_{k=0}^\infty \frac{-(x/4)^4}{[(3/2)_k]^2[(1)_k]^2}
   * @f]
   */
  template<typename _Real>
    _Real
    __kelvin_bex_series(_Real __x, int __sign)
    {
      using _WijnSum = __gnu_cxx::_VanWijngaardenSum<_Real>;
      using _WenigerDeltaWijnSum = __gnu_cxx::_WenigerDeltaSum<_WijnSum>;

      constexpr auto _S_1d2 = _Real{1} / _Real{2};
      constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __xd4 = __x / _Real{4};
      const auto __tmp = __xd4 * __xd4;
      const auto __y = __tmp * __tmp;
      auto __term = _Real{1};
      _WenigerDeltaWijnSum __bex;
      __bex += __term;
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  auto __fact = __k * (__k + __sign * _S_1d2);
	  __term *= -__y / __fact / __fact;
	  __bex += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__bex()))
	    break;
	}
      return __bex();
    }


  /**
   *
   */
  template<typename _Real>
    _Real
    __kelvin_kex_series(_Real __x, int __sign)
    {
      using _BasicSum = __gnu_cxx::_BasicSum<_Real>;
      using _AitkenSum = __gnu_cxx::_AitkenDeltaSquaredSum<_BasicSum>;
      using _WijnSum = __gnu_cxx::_VanWijngaardenSum<_Real>;
      using _WenigerDeltaWijnSum = __gnu_cxx::_WenigerDeltaSum<_WijnSum>;

      constexpr auto _S_1d2 = _Real{1} / _Real{2};
      constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __xd4 = __x / _Real{4};
      const auto __tmp = __xd4 * __xd4;
      const auto __y = __tmp * __tmp;
      auto __term = _Real{1};
      _AitkenSum _H_n;
      _H_n += 1;
      _WenigerDeltaWijnSum __kex;
      __kex += __term;
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  auto __fact = __k * (__k + __sign * _S_1d2);
	  __term *= -__y / __fact / __fact;

	  _H_n += _Real{1} / _Real(2 * __k);
	  auto __hterm = _Real{1} / _Real(2 * __k + 1);
	  if (__sign == +1)
	    _H_n += __hterm;
	  __kex += __term * _H_n();
	  if (__sign == -1)
	    _H_n += __hterm;
	  if (std::abs(__term) < _S_eps * std::abs(__kex()))
	    break;
	}
      return __kex();
    }


  /**
   * Return the Kelvin function @f$ ber(x) @f$ for real argument @c x.
   */
  template<typename _Real>
    inline _Real
    __kelvin_ber_series(_Real __x)
    { return __kelvin_bex_series(__x, -1); }


  /**
   * Return the Kelvin function @f$ bei(x) @f$ for real argument @c x.
   */
  template<typename _Real>
    inline _Real
    __kelvin_bei_series(_Real __x)
    { return __x * __x * __kelvin_bex_series(__x, +1) / _Real{4}; }


  /**
   * Return the irregular Kelvin function @f$ ker(x) @f$ for real argument @c x.
   */
  template<typename _Real>
    _Real
    __kelvin_ker_series(_Real __x)
    {
      constexpr auto _S_gamma_e = __gnu_cxx::__math_constants<_Real>::__gamma_e;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Real>::__pi_quarter;
      auto __ker = __kelvin_kex_series(__x, -1);
      auto __ber = __kelvin_bex_series(__x, -1);
      auto __bei = __kelvin_bex_series(__x, +1);
      auto __xxd4 = __x * __x / _Real{4};
      auto __ln = std::log(__x / _Real{2}) + _S_gamma_e;
      return -__ln * __ber + _S_pi_4 * __xxd4 * __bei + __ker - _Real{1};
    }


  /**
   * Return the irregular Kelvin function @f$ ker(x) @f$ for real argument @c x.
   */
  template<typename _Real>
    _Real
    __kelvin_kei_series(_Real __x)
    {
      constexpr auto _S_gamma_e = __gnu_cxx::__math_constants<_Real>::__gamma_e;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Real>::__pi_quarter;
      auto __kei = __kelvin_kex_series(__x, +1);
      auto __ber = __kelvin_bex_series(__x, -1);
      auto __bei = __kelvin_bex_series(__x, +1);
      auto __xxd4 = __x * __x / _Real{4};
      auto __ln = std::log(__x / _Real{2}) + _S_gamma_e;
      return -__ln * __xxd4 * __bei - _S_pi_4 * __ber + __xxd4 * __kei;
    }


  /**
   * Compute the Kelvin functions by series expansion.
   * The Kelvin functions and their expansions are:
   * @f[
   *    ber(x) = \sum_{k=0}^\infty  \frac{(-x^4/16)}{[(2k)!]^2}
   * @f]
   * @f[
   *    bei(x) = \frac{x^2}{4}\sum_{k=0}^\infty \frac{(-x^4/16)}{[(2k+1)!]^2}
   * @f]
   * @f[
   *    ker(x) = -\left[log\left(\frac{x}{2}\right) + \gamma_E\right] ber(x)
   *           + \frac{\pi}{4} bei(x)
   *           + \sum_{k=0}^\infty \frac{(-x^4/16)}{[(2k)!]^2}H_{2k}
   * @f]
   * @f[
   *    kei(x) = -\left[log\left(\frac{x}{2}\right) + \gamma_E\right] bei(x)
   *           - \frac{\pi}{4} ber(x)
   *           + \frac{x^2}{4}
   *               \sum_{k=0}^\infty \frac{(-x^4/16)}{[(2k+1)!]^2}H_{2k+1}
   * @f]
   * where
   * @f[
   *    H_n = \sum_{k=1}^n \frac{1}{k}
   * @f]
   * is the harmonic number.
   */
  template<typename _Real>
    _KelvinState<_Real>
    __kelvin_series(_Real __x)
    {
      using _BasicSum = __gnu_cxx::_BasicSum<_Real>;
      using _AitkenSum = __gnu_cxx::_AitkenDeltaSquaredSum<_BasicSum>;
      using _WijnSum = __gnu_cxx::_VanWijngaardenSum<_Real>;
      using _WenigerDeltaWijnSum = __gnu_cxx::_WenigerDeltaSum<_WijnSum>;

      constexpr auto _S_gamma_e = __gnu_cxx::__math_constants<_Real>::__gamma_e;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Real>::__pi_quarter;
      constexpr auto _S_eps = _Real{0.01L}
				* std::numeric_limits<_Real>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __xd2 = __x / _Real{2};
      const auto __xxd4 = __xd2 * __xd2;
      const auto __y = __xxd4 * __xxd4;
      auto __termr = _Real{1};
      auto __termi = _Real{1};
      _AitkenSum _H_n;
      _WenigerDeltaWijnSum __ber, __bei, __ker, __kei;
      __ber += __termr;
      __ker += __termr;
      __bei += __termi;
      __kei += __termi;
      _H_n += _Real{1};
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  auto __factr = _Real{1} / _Real(2 * __k - 1) / _Real(2 * __k);
	  __termr *= -__y * __factr * __factr;
	  __ber += __termr;

	  _H_n += _Real{1} / _Real(2 * __k);
	  __ker += __termr * _H_n();

	  auto __facti = _Real{1} / _Real(2 * __k) / _Real(2 * __k + 1);
	  __termi *= -__y * __facti * __facti;
	  __bei += __termi;

	  _H_n += _Real{1} / _Real(2 * __k + 1);
	  __kei += __termi * _H_n();

	  if (std::abs(__termr) < _S_eps * std::abs(__ber()))
	    break;
	}
      auto __ln = std::log(__x / _Real{2}) + _S_gamma_e;
      return _KelvinState<_Real>{_Real{0}, __x, __ber(), __xxd4 * __bei(),
	      -__ln * __ber() + _S_pi_4 * __xxd4 * __bei() + __ker() - _Real{1},
	      -__ln * __xxd4 * __bei() - _S_pi_4 * __ber() + __xxd4 * __kei()};
    }


  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename _Real>
    _KelvinState<_Real>
    __kelvin_asymp(_Real __x)
    {
      using _Cmplx = std::complex<_Real>;
      using _BasicSum = __gnu_cxx::_BasicSum<_Cmplx>;
      using _WenigerDeltaSum = __gnu_cxx::_WenigerDeltaSum<_BasicSum>;

      constexpr auto _S_j = _Cmplx{0, 1};
      constexpr auto _S_1d2 = _Real{1} / _Real{2};
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Real>::__pi_quarter;
      constexpr auto _S_pi_8 = _S_1d2 * _S_pi_4;
      constexpr auto _S_3pi_4 = _Real{3} * _S_pi_4;
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Real>::__root_2;
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Real>::__root_pi;
      constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __y = _Real{1} / (_Real{32} * __x);
      auto __term = _Real{1};
      const auto __xrt2 = __x / _S_sqrt_2;
      auto __barg = __xrt2 - _S_pi_8;
      auto __karg = -__xrt2 - _S_pi_8;
      _WenigerDeltaSum __be, __ke;
      __be += std::polar(__term, __barg);
      __ke += std::polar(__term, __karg);
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  __barg -= _S_pi_4;
	  __karg += _S_3pi_4;
	  auto __fact = _Real(4 * __k - 2);
	  auto __next = __y * __fact * __fact / _Real(__k);
	  if (std::abs(__next) > _Real{1})
	    break;
	  __term *= __next;
	  __be += std::polar(__term, __barg);
	  __ke += std::polar(__term, __karg);

	  if (std::abs(__term) < _S_eps * std::abs(__be()))
	    break;
	}
      const auto __exp = std::exp(__xrt2);
      const auto __rt = std::sqrt(_Real{2} * __x);
      const auto __kfact = _S_sqrt_pi / __rt / __exp;
      const auto __kex = __kfact * __ke();
      const auto __bfact = __exp / _S_sqrt_pi / __rt;
      const auto __bex = __bfact * __be() + _S_j * __kex / _S_pi;
      return _KelvinState<_Real>{_Real{0}, __x,
			         std::real(__bex), std::imag(__bex),
			         std::real(__kex), std::imag(__kex)};
    }

  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename _Real>
    inline _Real
    __kelvin_ber_asymp(_Real __x)
    { return __kelvin_asymp(__x).__ber; }

  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename _Real>
    inline _Real
    __kelvin_bei_asymp(_Real __x)
    { return __kelvin_asymp(__x).__bei; }

  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename _Real>
    inline _Real
    __kelvin_ker_asymp(_Real __x)
    { return __kelvin_asymp(__x).__ker; }

  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename _Real>
    inline _Real
    __kelvin_kei_asymp(_Real __x)
    { return __kelvin_asymp(__x).__kei; }


  /**
   * Compute the Kelvin functions of integer order n by series expansion.
   */
  template<typename _Real>
    _KelvinState<_Real>
    __kelvin_series(int __n, _Real __x)
    {
      using _Cmplx = std::complex<_Real>;
      constexpr auto _S_j = _Cmplx{0, 1};
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Real>::__pi_half;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Real>::__pi_quarter;
      constexpr auto _S_3pi_4 = _Real{3} * _S_pi_4;
      constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();
      constexpr auto _S_maxiter = 1000;
      if (__n < 0)
	{
	  const auto _Cnp = _Real(1 - 2 * (__n & 1));
	  auto _Kv = __kelvin_series(-__n, __x);
	  return _KelvinState<_Real>{_Real(__n), __x,
	  			     _Cnp * _Kv.__ber, _Cnp * _Kv.__bei,
				     _Cnp * _Kv.__ker, _Cnp * _Kv.__kei};
	}
      else
	{
	  const auto __xd2 = __x / _Real{2};
	  const auto __y = __xd2 * __xd2;

	  __gnu_cxx::_BasicSum<_Cmplx> __be;
	  auto __bterm = _Real{1};
	  auto __barg = _S_3pi_4 * _Real(__n);
	  __be += std::polar(__bterm, __barg);
	  for (auto __k = 1; __k < _S_maxiter; ++__k)
	    {
	      __bterm *= __y / _Real(__k + __n) / _Real(__k);
	      __barg += _S_pi_2;
	      __be += std::polar(__bterm, __barg);
	      if (std::abs(__bterm) < _S_eps * std::abs(__be()))
		break;
	    }

	  __gnu_cxx::_BasicSum<_Cmplx> __ke1;
	  if (__n > 0)
	    {
	      auto __kterm1 = _Real{1};
	      auto __karg1 = _S_3pi_4 * _Real(__n);
	      __ke1 += std::polar(__kterm1, -__karg1);
	      for (auto __k = 1; __k < __n - 1; ++__k)
		{
		  __kterm1 *= __y / _Real(__n - 1 - __k) / _Real(__k);
		  __karg1 += _S_pi_2;
		  __ke1 += std::polar(__kterm1, -__karg1);
		}
	      if (__n > 1)
		{
		  __kterm1 *= __y / _Real(__n - 1);
		  __karg1 += _S_pi_2;
		  __ke1 += std::polar(__kterm1, -__karg1);
		}
	    }

	  __gnu_cxx::_BasicSum<_Cmplx> __ke2;
	  auto __hsum2 = __psi<_Real>(1) + __psi<_Real>(1 + __n);
	  auto __kterm2 = _Real{1};
	  auto __karg2 = _S_3pi_4 * _Real(__n);
	  __ke2 += std::polar(__hsum2 * __kterm2, __karg2);
	  for (auto __k = 1; __k < _S_maxiter; ++__k)
	    {
	      __hsum2 += _Real{1} / _Real(1 + __k)
		       + _Real{1} / _Real(1 + __n + __k);
	      __kterm2 *= __y / _Real(__k) / _Real(__n + __k);
	      __karg2 += _S_pi_2;
	      __ke2 += std::polar(__hsum2 * __kterm2, __karg2);
	      if (std::abs(__hsum2 * __kterm2) < _S_eps * std::abs(__ke2()))
		break;
	    }

	  const auto __nfact = __factorial<_Real>(__n);
	  const auto __pow = std::pow(__xd2, _Real(__n));
	  auto __bex = __pow * __be() / __nfact;
	  auto __kex = -std::log(__xd2) * __bex - _S_j * _S_pi_4 * __bex
		     + __ke1() / __pow / _Real{2}
		       / (__n > 0 ? __factorial<_Real>(__n - 1) : _Real{1})
		     + __pow * __ke2() / __nfact / _Real{2}
		     - _Real{1}; // I needed this for ker(x) too.
	  return _KelvinState<_Real>{_Real(__n), __x,
				     std::real(__bex), std::imag(__bex),
				     std::real(__kex), std::imag(__kex)};
	}
    }


  /**
   * Compute the Kelvin functions of order @f$ \nu @f$ by series expansion.
   */
  template<typename _Real>
    _KelvinState<_Real>
    __kelvin_series(_Real __nu, _Real __x)
    {
      using _Cmplx = std::complex<_Real>;
      constexpr auto _S_j = _Cmplx{0, 1};
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Real>::__pi_half;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Real>::__pi_quarter;
      constexpr auto _S_3pi_4 = _Real{3} * _S_pi_4;
      constexpr auto _S_maxiter = 1000;

      // I have to try C++17 init-select!
      if (auto __isint = __gnu_cxx::__fpinteger(__nu); __isint)
	  return __kelvin_series(__isint.__value, __x);
/*
      if (__gnu_cxx::__fpinteger(__nu))
	{
	  int __n = std::nearbyint(__nu);
	  return __kelvin_series(__n, __x);
	}
*/
      else
	{
	  const auto __xd2 = __x / _Real{2};
	  const auto __y = __xd2 * __xd2;

	  __gnu_cxx::_BasicSum<_Cmplx> __bep;
	  auto __termp = _Real{1};
	  auto __argp = _S_3pi_4 * __nu;
	  __bep += std::polar(__termp, __argp);
	  for (auto __k = 1; __k < _S_maxiter; ++__k)
	    {
	      __termp *= __y / _Real(__k + __nu) / _Real(__k);
	      __argp += _S_pi_2;
	      __bep += std::polar(__termp, __argp);
	    }
	  auto __bex = std::pow(__xd2, __nu) * __bep()
		     / std::tgamma(_Real{1} + __nu);

	  __gnu_cxx::_BasicSum<_Cmplx> __bem;
	  auto __termm = _Real{1};
	  auto __argm = -_S_3pi_4 * __nu;
	  __bem += std::polar(__termm, __argm);
	  for (auto __k = 1; __k < _S_maxiter; ++__k)
	    {
	      __termm *= __y / _Real(__k - __nu) / _Real(__k);
	      __argm += _S_pi_2;
	      __bem += std::polar(__termm, __argm);
	    }
	  auto __bey = std::pow(__xd2, -__nu) * __bem()
		     / std::tgamma(_Real{1} - __nu);
	  const auto __nupi = __nu * _S_pi;
	  const auto __csc = _Real{1} / std::sin(__nupi);
	  const auto __cot = _Real{1} / std::tan(__nupi);
	  auto __kex = _S_pi_2 * (__csc * __bey - __cot * __bex + _S_j * __bex);

	  return _KelvinState<_Real>{__nu, __x,
				     std::real(__bex), std::imag(__bex),
				     std::real(__kex), std::imag(__kex)};
	}
    }


  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename _Real>
    _KelvinState<_Real>
    __kelvin_asymp(_Real __nu, _Real __x)
    {
      using _Cmplx = std::complex<_Real>;
      using _BasicSum = __gnu_cxx::_BasicSum<_Cmplx>;
      using _WenigerDeltaSum = __gnu_cxx::_WenigerDeltaSum<_BasicSum>;

      constexpr auto _S_j = _Cmplx{0, 1};
      constexpr auto _S_1d2 = _Real{1} / _Real{2};
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Real>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Real>::__pi_half;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Real>::__pi_quarter;
      constexpr auto _S_pi_8 = _S_1d2 * _S_pi_4;
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Real>::__root_2;
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Real>::__root_pi;
      constexpr auto _S_eps = std::numeric_limits<_Real>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __y = _Real{1} / (_Real{2} * __x);
      auto __bterm = _Real{1};
      auto __kterm = _Real{1};
      const auto __xrt2 = __x / _S_sqrt_2;
      auto __barg = __xrt2 + __nu * _S_pi_2 - _S_pi_8;
      auto __karg = __xrt2 + __nu * _S_pi_2 + _S_pi_8;
      _WenigerDeltaSum __be, __ke;
      __be += std::polar(__bterm, __barg);
      __ke += std::polar(__kterm, __karg);
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  __barg -= _S_pi_4;
	  __karg += _S_pi_4;
	  auto __fact = (_Real(__k) - _Real{0.5L} - __nu)
		      * (_Real(__k) - _Real{0.5L} + __nu) / _Real(__k);
	  auto __next = __y * __fact;
	  if (std::abs(__next) > _Real{1})
	    break;
	  __bterm *= __next;
	  __kterm *= -__bterm;
	  __be += std::polar(__bterm, __barg);
	  __ke += std::polar(__kterm, __karg);

	  if (std::abs(__bterm) < _S_eps * std::abs(__be()))
	    break;
	}
      const auto __exp = std::exp(__xrt2);
      const auto __rt = _S_sqrt_2 * std::sqrt(__x);
      const auto __kfact = _S_sqrt_pi / __rt / __exp;
      const auto __kex = __kfact * __ke();
      const auto __bfact = __exp / _S_sqrt_pi / __rt;
      const auto __bex = __bfact * __be() + _S_j * __kex / _S_pi;
      return _KelvinState<_Real>{__nu, __x,
			         std::real(__bex), std::imag(__bex),
			         std::real(__kex), -std::imag(__kex)};
    }

  /**
   * Compute the Kelvin function @f$ ber(x)@f$.
   */
  template<typename _Real>
    inline _Real
    __kelvin_ber(_Real __x)
    {
      const auto _S_switch = _Real{26};
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Real>();
      else if (std::abs(__x) < _S_switch)
	return __kelvin_ber_series(__x);
      else
	return __kelvin_ber_asymp(__x);
    }

  /**
   * Compute the Kelvin function @f$ bei(x)@f$.
   */
  template<typename _Real>
    inline _Real
    __kelvin_bei(_Real __x)
    {
      const auto _S_switch = _Real{26};
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Real>();
      else if (std::abs(__x) < _S_switch)
	return __kelvin_bei_series(__x);
      else
	return __kelvin_bei_asymp(__x);
    }

  /**
   * Compute the Kelvin function @f$ ker(x)@f$.
   */
  template<typename _Real>
    inline _Real
    __kelvin_ker(_Real __x)
    {
      const auto _S_switch = _Real{5};
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Real>();
      else if (std::abs(__x) < _S_switch)
	return __kelvin_ker_series(__x);
      else
	return __kelvin_ker_asymp(__x);
    }

  /**
   * Compute the Kelvin function @f$ kei(x)@f$.
   */
  template<typename _Real>
    inline _Real
    __kelvin_kei(_Real __x)
    {
      const auto _S_switch = _Real{5};
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Real>();
      else if (std::abs(__x) < _S_switch)
	return __kelvin_kei_series(__x);
      else
	return __kelvin_kei_asymp(__x);
    }

  /**
   * Compute the Kelvin function @f$ ber_\nu(x)@f$.
   */
  template<typename _Real>
    inline _Real
    __kelvin_ber(_Real __nu, _Real __x)
    {
      const auto _S_switch = _Real{26};
      if (__isnan(__nu) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Real>();
      else if (std::abs(__x) < _S_switch)
	return __kelvin_series(__nu, __x).__ber;
      else
	return __kelvin_asymp(__nu, __x).__ber;
    }

  /**
   * Compute the Kelvin function @f$ bei_\nu(x)@f$.
   */
  template<typename _Real>
    inline _Real
    __kelvin_bei(_Real __nu, _Real __x)
    {
      const auto _S_switch = _Real{26};
      if (__isnan(__nu) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Real>();
      else if (std::abs(__x) < _S_switch)
	return __kelvin_series(__nu, __x).__bei;
      else
	return __kelvin_asymp(__nu, __x).__bei;
    }

  /**
   * Compute the Kelvin function @f$ ker_\nu(x)@f$.
   */
  template<typename _Real>
    inline _Real
    __kelvin_ker(_Real __nu, _Real __x)
    {
      const auto _S_switch = _Real{5};
      if (__isnan(__nu) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Real>();
      else if (std::abs(__x) < _S_switch)
	return __kelvin_series(__nu, __x).__ker;
      else
	return __kelvin_asymp(__nu, __x).__ker;
    }

  /**
   * Compute the Kelvin function @f$ kei_\nu(x)@f$.
   */
  template<typename _Real>
    inline _Real
    __kelvin_kei(_Real __nu, _Real __x)
    {
      const auto _S_switch = _Real{5};
      if (__isnan(__nu) || __isnan(__x))
	return __gnu_cxx::__quiet_NaN<_Real>();
      else if (std::abs(__x) < _S_switch)
	return __kelvin_series(__nu, __x).__kei;
      else
	return __kelvin_asymp(__nu, __x).__kei;
    }

} // namespace __detail
} // namespace std

namespace __gnu_cxx
{


  //  Kelvin functions

  /**
   * Return the Kelvin function @f$ ber(x) @f$ for @c float argument @c x.
   *
   * @see kelvin_ber for details.
   */
  inline float
  kelvin_berf(float __x)
  { return std::__detail::__kelvin_ber<float>(__x); }

  /**
   * Return the Kelvin function @f$ ber(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see kelvin_ber for details.
   */
  inline long double
  kelvin_berl(long double __x)
  { return std::__detail::__kelvin_ber<long double>(__x); }

  /**
   * Return the Kelvin function @f$ ber(x) @f$ for @c real argument @c x.
   *
   * Kelvin integrals @f$ ber(x) @f$ and @f$ bei(x) @f$ are given by
   * \f[
   *   J_0(\frac{x-ix}{\sqrt{2}}) = ber(x) + i bei(x)
   * \f]
   * where @f$ J_0(x) @f$ is the Bessel function of the first kind.
   *
   * @tparam _Real The floating-point type of the argument @c __x.
   * @param  __x  The argument of the Kelvin function
   */
  template<typename _Real>
    inline typename __gnu_cxx::__promote<_Real>::__type
    kelvin_ber(_Real __x)
    {
      typedef typename __gnu_cxx::__promote<_Real>::__type __type;
      return std::__detail::__kelvin_ber<__type>(__x);
    }

  /**
   * Return the Kelvin function @f$ bei(x) @f$ for @c float argument @c x.
   *
   * @see kelvin_bei for details.
   */
  inline float
  kelvin_beif(float __x)
  { return std::__detail::__kelvin_bei<float>(__x); }

  /**
   * Return the Kelvin function @f$ bei(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see kelvin_bei for details.
   */
  inline long double
  kelvin_beil(long double __x)
  { return std::__detail::__kelvin_bei<long double>(__x); }

  /**
   * Return the Kelvin function @f$ bei(x) @f$ for @c real argument @c x.
   *
   * Kelvin integrals @f$ bei(x) @f$ and @f$ bei(x) @f$ are given by
   * \f[
   *   J_0(\frac{x-ix}{\sqrt{2}}) = ber(x) + i bei(x)
   * \f]
   * where @f$ J_0(x) @f$ is the Bessel function of the first kind.
   *
   * @tparam _Real The floating-point type of the argument @c __x.
   * @param  __x  The argument of the Kelvin function
   */
  template<typename _Real>
    inline typename __gnu_cxx::__promote<_Real>::__type
    kelvin_bei(_Real __x)
    {
      typedef typename __gnu_cxx::__promote<_Real>::__type __type;
      return std::__detail::__kelvin_bei<__type>(__x);
    }

  /**
   * Return the Kelvin function @f$ ker(x) @f$ for @c float argument @c x.
   *
   * @see kelvin_ker for details.
   */
  inline float
  kelvin_kerf(float __x)
  { return std::__detail::__kelvin_ker<float>(__x); }

  /**
   * Return the Kelvin function @f$ ker(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see kelvin_ker for details.
   */
  inline long double
  kelvin_kerl(long double __x)
  { return std::__detail::__kelvin_ker<long double>(__x); }

  /**
   * Return the Kelvin function @f$ ker(x) @f$ for @c real argument @c x.
   *
   * Kelvin integrals @f$ ker(x) @f$ and @f$ kei(x) @f$ are given by
   * \f[
   *   J_0(\frac{x-ix}{\sqrt{2}}) = ker(x) + i kei(x)
   * \f]
   * where @f$ J_0(x) @f$ is the Bessel function of the first kind.
   *
   * @tparam _Real The floating-point type of the argument @c __x.
   * @param  __x  The argument of the Kelvin function
   */
  template<typename _Real>
    inline typename __gnu_cxx::__promote<_Real>::__type
    kelvin_ker(_Real __x)
    {
      typedef typename __gnu_cxx::__promote<_Real>::__type __type;
      return std::__detail::__kelvin_ker<__type>(__x);
    }

  /**
   * Return the Kelvin function @f$ kei(x) @f$ for @c float argument @c x.
   *
   * @see kelvin_kei for details.
   */
  inline float
  kelvin_keif(float __x)
  { return std::__detail::__kelvin_kei<float>(__x); }

  /**
   * Return the Kelvin function @f$ kei(x) @f$
   * for <tt>long double</tt> argument @c x.
   *
   * @see kelvin_kei for details.
   */
  inline long double
  kelvin_keil(long double __x)
  { return std::__detail::__kelvin_kei<long double>(__x); }

  /**
   * Return the Kelvin function @f$ kei(x) @f$ for @c real argument @c x.
   *
   * Kelvin integrals @f$ kei(x) @f$ and @f$ kei(x) @f$ are given by
   * \f[
   *   J_0(\frac{x-ix}{\sqrt{2}}) = ker(x) + i kei(x)
   * \f]
   * where @f$ J_0(x) @f$ is the Bessel function of the first kind.
   *
   * @tparam _Real The floating-point type of the argument @c __x.
   * @param  __x  The argument of the Kelvin function
   */
  template<typename _Real>
    inline typename __gnu_cxx::__promote<_Real>::__type
    kelvin_kei(_Real __x)
    {
      typedef typename __gnu_cxx::__promote<_Real>::__type __type;
      return std::__detail::__kelvin_kei<__type>(__x);
    }

  /**
   * Return the Kelvin function @f$ ber_\nu(x) @f$ for float
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_ber for details.
   */
  inline float
  kelvin_berf(float __nu, float __x)
  { return std::__detail::__kelvin_ber<float>(__nu, __x); }

  /**
   * Return the Kelvin function @f$ ber_\nu(x) @f$ for <tt>long double</tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_ber for details.
   */
  inline long double
  kelvin_berl(long double __nu, long double __x)
  { return std::__detail::__kelvin_ber<long double>(__nu, __x); }

  /**
   * Return the Kelvin function @f$ ber_\nu(x) @f$ of real
   * order @f$ \nu @f$ and argument @f$ x @f$.
   * The Kelvin functions @f$ ber_\nu(x) @f$ and @f$ bei_\nu(x) @f$
   * are defined by:
   * @f[
   *    J_\nu(\frac{x-ix}{\sqrt{2}}) = J_\nu(xe^{i3\pi/4})
   *     = ber_\nu(x) + i bei_\nu(x)
   * @f]
   * where J_\nu(x) is the Bessel function of the first kind.
   *
   * @tparam _Real The real type of the argument
   * @param __nu The real order
   * @param __x The real argument
   */
  template<typename _Real>
    inline __gnu_cxx::__promote_fp_t<_Real>
    kelvin_ber(_Real __nu, _Real __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Real>;
      return std::__detail::__kelvin_ber<__type>(__nu, __x);
    }

  /**
   * Return the Kelvin function @f$ bei_\nu(x) @f$ for float
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_bei for details.
   */
  inline float
  kelvin_beif(float __nu, float __x)
  { return std::__detail::__kelvin_bei<float>(__nu, __x); }

  /**
   * Return the Kelvin function @f$ bei_\nu(x) @f$ for <tt>long double</tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_bei for details.
   */
  inline long double
  kelvin_beil(long double __nu, long double __x)
  { return std::__detail::__kelvin_bei<long double>(__nu, __x); }

  /**
   * Return the Kelvin function @f$ bei_\nu(x) @f$ of real
   * order @f$ \nu @f$ and argument @f$ x @f$.
   * The Kelvin functions @f$ bei_\nu(x) @f$ and @f$ bei_\nu(x) @f$
   * are defined by:
   * @f[
   *    J_\nu(\frac{x-ix}{\sqrt{2}}) = J_\nu(xe^{i3\pi/4})
   *     = ber_\nu(x) + i bei_\nu(x)
   * @f]
   * where J_\nu(x) is the Bessel function of the first kind.
   *
   * @tparam _Real The real type of the argument
   * @param __nu The real order
   * @param __x The real argument
   */
  template<typename _Real>
    inline __gnu_cxx::__promote_fp_t<_Real>
    kelvin_bei(_Real __nu, _Real __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Real>;
      return std::__detail::__kelvin_bei<__type>(__nu, __x);
    }

  /**
   * Return the Kelvin function @f$ ker_\nu(x) @f$ for float
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_ker for details.
   */
  inline float
  kelvin_kerf(float __nu, float __x)
  { return std::__detail::__kelvin_ker<float>(__nu, __x); }

  /**
   * Return the Kelvin function @f$ ker_\nu(x) @f$ for <tt>long double</tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_ker for details.
   */
  inline long double
  kelvin_kerl(long double __nu, long double __x)
  { return std::__detail::__kelvin_ker<long double>(__nu, __x); }

  /**
   * Return the Kelvin function @f$ ker_\nu(x) @f$ of real
   * order @f$ \nu @f$ and argument @f$ x @f$.
   * The Kelvin functions @f$ ker_\nu(x) @f$ and @f$ kei_\nu(x) @f$
   * are defined by:
   * @f[
   *    K_\nu(\frac{x+ix}{\sqrt{2}}) = K_\nu(xe^{i\pi/4})
   *     = ker_\nu(x) + i kei_\nu(x)
   * @f]
   * where K_\nu(x) is the regular modified Bessel function.
   *
   * @tparam _Real The real type of the argument
   * @param __nu The real order
   * @param __x The real argument
   */
  template<typename _Real>
    inline __gnu_cxx::__promote_fp_t<_Real>
    kelvin_ker(_Real __nu, _Real __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Real>;
      return std::__detail::__kelvin_ker<__type>(__nu, __x);
    }

  /**
   * Return the Kelvin function @f$ kei_\nu(x) @f$ for float
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_kei for details.
   */
  inline float
  kelvin_keif(float __nu, float __x)
  { return std::__detail::__kelvin_kei<float>(__nu, __x); }

  /**
   * Return the Kelvin function @f$ kei_\nu(x) @f$ for <tt>long double</tt>
   * order @f$ \nu @f$ and argument @c x.
   *
   * @see kelvin_kei for details.
   */
  inline long double
  kelvin_keil(long double __nu, long double __x)
  { return std::__detail::__kelvin_kei<long double>(__nu, __x); }

  /**
   * Return the Kelvin function @f$ kei_\nu(x) @f$ of real
   * order @f$ \nu @f$ and argument @f$ x @f$.
   * The Kelvin functions @f$ kei_\nu(x) @f$ and @f$ kei_\nu(x) @f$
   * are defined by:
   * @f[
   *    K_\nu(\frac{x+ix}{\sqrt{2}}) = K_\nu(xe^{i\pi/4})
   *     = ker_\nu(x) + i kei_\nu(x)
   * @f]
   * where K_\nu(x) is the regular modified Bessel function.
   *
   * @tparam _Real The real type of the argument
   * @param __nu The real order
   * @param __x The real argument
   */
  template<typename _Real>
    inline __gnu_cxx::__promote_fp_t<_Real>
    kelvin_kei(_Real __nu, _Real __x)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Real>;
      return std::__detail::__kelvin_kei<__type>(__nu, __x);
    }

} // namespace __gnu_cxx


/**
 * Run Kelvin with individual sums: __kelvin_ber_series, etc.
 */
template<typename _Real>
  void
  run_kelvin1()
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::__detail::__kelvin_ber_series(_Real{});

    std::cout << "\n\nPrint Kelvin functions computed by series expansions\n";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    std::cout << "\n\n";
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "ber"
	      << std::setw(width) << "bei"
	      << std::setw(width) << "ker"
	      << std::setw(width) << "kei"
	      << '\n';
    std::cout << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << '\n';
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _Real(0.1L) * i;
	auto ber = std::__detail::__kelvin_ber_series(x);
	auto bei = std::__detail::__kelvin_bei_series(x);
	auto ker = std::__detail::__kelvin_ker_series(x);
	auto kei = std::__detail::__kelvin_kei_series(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << ber
		  << std::setw(width) << bei
		  << std::setw(width) << ker
		  << std::setw(width) << kei
		  << '\n';
      }
    std::cout << std::endl;
  }


/**
 * Run Kelvin with __kelvin_series which combines sums.
 */
template<typename _Real>
  void
  run_kelvin2()
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::__detail::__kelvin_series(_Real{});

    std::cout << "\n\nPrint Kelvin functions computed by series expansions\n";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    std::cout << "\n\n";
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "ber"
	      << std::setw(width) << "bei"
	      << std::setw(width) << "ker"
	      << std::setw(width) << "kei"
	      << '\n';
    std::cout << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << '\n';
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _Real(0.1L) * i;
	auto ke = std::__detail::__kelvin_series(x);
	std::cout << std::setw(width) << ke.__x
		  << std::setw(width) << ke.__ber
		  << std::setw(width) << ke.__bei
		  << std::setw(width) << ke.__ker
		  << std::setw(width) << ke.__kei
		  << '\n';
      }
    std::cout << std::endl;
  }


/**
 * Diff Kelvin __kelvin_asymp with __kelvin_series.
 */
template<typename _Real>
  void
  diff_kelvin2()
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "\n\nDiff Kelvin functions computed by series expansions\n";
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    std::cout << "\n\n";
    std::cout << "  " << std::setw(width) << "x"
	      << "  " << std::setw(width) << "ber asymp"
	      << "  " << std::setw(width) << "ber series"
	      << "  " << std::setw(width) << "bei asymp"
	      << "  " << std::setw(width) << "bei series"
	      << "  " << std::setw(width) << "ker asymp"
	      << "  " << std::setw(width) << "ker series"
	      << "  " << std::setw(width) << "kei asymp"
	      << "  " << std::setw(width) << "kei series"
	      << '\n';
    std::cout << "  " << std::setw(width) << "============="
	      << "  " << std::setw(width) << "============="
	      << "  " << std::setw(width) << "============="
	      << "  " << std::setw(width) << "============="
	      << "  " << std::setw(width) << "============="
	      << "  " << std::setw(width) << "============="
	      << "  " << std::setw(width) << "============="
	      << "  " << std::setw(width) << "============="
	      << "  " << std::setw(width) << "============="
	      << '\n';
    for (int i = 50; i <= 400; ++i)
      {
	auto x = _Real(0.1L) * i;
	auto kes = std::__detail::__kelvin_series(x);
	auto kea = std::__detail::__kelvin_asymp(x);
	std::cout << "  " << std::setw(width) << kes.__x
		  << "  " << std::setw(width) << kea.__ber
		  << "  " << std::setw(width) << kes.__ber
		  << "  " << std::setw(width) << kea.__bei
		  << "  " << std::setw(width) << kes.__bei
		  << "  " << std::setw(width) << kea.__ker
		  << "  " << std::setw(width) << kes.__ker
		  << "  " << std::setw(width) << kea.__kei
		  << "  " << std::setw(width) << kes.__kei
		  << '\n';
      }
    std::cout << std::endl;
  }


/**
 * Run Kelvin with __kelvin_series(nu, x) which combines sums for general order.
 */
template<typename _Real>
  void
  run_kelvin3(_Real nu = _Real{0})
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "\n\nPrint Kelvin functions computed by series expansions, nu = " << nu << '\n';
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    std::cout << "\n\n";
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "ber"
	      << std::setw(width) << "bei"
	      << std::setw(width) << "ker"
	      << std::setw(width) << "kei"
	      << '\n';
    std::cout << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << '\n';
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _Real(0.1L) * i;
	auto ke = std::__detail::__kelvin_series(nu, x);
	std::cout << std::setw(width) << ke.__x
		  << std::setw(width) << ke.__ber
		  << std::setw(width) << ke.__bei
		  << std::setw(width) << ke.__ker
		  << std::setw(width) << ke.__kei
		  << '\n';
      }
    std::cout << std::endl;
  }


/**
 * Diff Kelvin __kelvin_asymp(nu, x) with __kelvin_series(nu, x) for general order.
 */
template<typename _Real>
  void
  diff_kelvin3(_Real nu = _Real{0})
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "\n\nDiff Kelvin functions computed by series expansions, nu = " << nu << '\n';
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    std::cout << "\n\n";
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "ber"
	      << std::setw(width) << "bei"
	      << std::setw(width) << "ker"
	      << std::setw(width) << "kei"
	      << '\n';
    std::cout << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << '\n';
    for (int i = 50; i <= 400; ++i)
      {
	auto x = _Real(0.1L) * i;
	auto kes = std::__detail::__kelvin_series(nu, x);
	auto kea = std::__detail::__kelvin_asymp(nu, x);
	std::cout << "  " << std::setw(width) << kes.__x
		  << "  " << std::setw(width) << (kea.__ber - kes.__ber) / std::abs(kes.__ber)
		  << "  " << std::setw(width) << (kea.__bei - kes.__bei) / std::abs(kes.__bei)
		  << "  " << std::setw(width) << (kea.__ker - kes.__ker) / std::abs(kes.__ker)
		  << "  " << std::setw(width) << (kea.__kei - kes.__kei) / std::abs(kes.__kei)
		  << '\n';
      }
    std::cout << std::endl;
  }


/**
 * Run Kelvin with __kelvin_series(n, x) which combines sums for integral order.
 */
template<typename _Real>
  void
  run_kelvin4(int n = 0)
  {
    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "\n\nPrint Kelvin functions computed by series expansions, n = " << n << '\n';
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";

    std::cout << "\n\n";
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "ber"
	      << std::setw(width) << "bei"
	      << std::setw(width) << "ker"
	      << std::setw(width) << "kei"
	      << '\n';
    std::cout << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << '\n';
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _Real(0.1L) * i;
	auto ke = std::__detail::__kelvin_series(n, x);
	std::cout << std::setw(width) << ke.__x
		  << std::setw(width) << ke.__ber
		  << std::setw(width) << ke.__bei
		  << std::setw(width) << ke.__ker
		  << std::setw(width) << ke.__kei
		  << '\n';
      }
    std::cout << std::endl;
  }


/**
 * Plot the scaled Kelvin functions.
 */
template<typename _Real>
  void
  plot_kelvin(std::string filename)
  {
    constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Real>::__root_2;
    constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Real>::__root_pi;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    data << "\n\n";
    data << "#"
	 << std::setw(width) << "x"
	 << std::setw(width) << "ber"
	 << std::setw(width) << "bei"
	 << std::setw(width) << "ker"
	 << std::setw(width) << "kei"
	 << std::setw(width) << "Wronskian"
	 << '\n';
    data << "#"
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << '\n';
    for (int i = 0; i <= +2000; ++i)
      {
	auto x = _Real(0.01Q * i);
	auto kv = std::__detail::__kelvin_series(x);
	auto exf = std::exp(x / _S_sqrt_2);
	auto rt2x = _S_sqrt_2 * std::sqrt(x);
	data << std::setw(width) << kv.__x
	     << std::setw(width) << _S_sqrt_pi * rt2x * kv.__ber / exf
	     << std::setw(width) << _S_sqrt_pi * rt2x * kv.__bei / exf
	     << std::setw(width) << rt2x * exf * kv.__ker / _S_sqrt_pi
	     << std::setw(width) << rt2x * exf * kv.__kei / _S_sqrt_pi
	     << '\n';
      }
    for (int i = 2001; i <= +4000; ++i)
      {
	auto x = _Real(0.01Q * i);
	auto kv = std::__detail::__kelvin_asymp(x);
	auto exf = std::exp(x / _S_sqrt_2);
	auto rt2x = _S_sqrt_2 * std::sqrt(x);
	data << std::setw(width) << kv.__x
	     << std::setw(width) << _S_sqrt_pi * rt2x * kv.__ber / exf
	     << std::setw(width) << _S_sqrt_pi * rt2x * kv.__bei / exf
	     << std::setw(width) << rt2x * exf * kv.__ker / _S_sqrt_pi
	     << std::setw(width) << rt2x * exf * kv.__kei / _S_sqrt_pi
	     << '\n';
      }
    data << "\n\n";
  }


/**
 * Plot the scaled Kelvin functions for orders 0, 1/2, 1, 3/2, 2, 5.
 */
template<typename _Real>
  void
  plot_kelvin_order(std::string filename)
  {
    constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Real>::__root_2;
    constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Real>::__root_pi;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    data << "\n\n";
    data << "#"
	 << std::setw(width) << "x"
	 << std::setw(width) << "ber"
	 << std::setw(width) << "bei"
	 << std::setw(width) << "ker"
	 << std::setw(width) << "kei"
	 << std::setw(width) << "Wronskian"
	 << '\n';
    data << "#"
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << '\n';

    std::vector<_Real> nuv{_Real{0.0L}, _Real{0.5L}, _Real{1.0L},
			   _Real{1.5L}, _Real{2.0L}, _Real{5.0L}};
    for (auto nu : nuv)
      {
	for (int i = 0; i <= 2000; ++i)
	  {
	    auto x = _Real(0.01Q * i);
	    auto kv = std::__detail::__kelvin_series(nu, x);
	    auto exf = std::exp(x / _S_sqrt_2);
	    auto rt2x = _S_sqrt_2 * std::sqrt(x);
	    data << std::setw(width) << kv.__x
		 << std::setw(width) << _S_sqrt_pi * rt2x * kv.__ber / exf
		 << std::setw(width) << _S_sqrt_pi * rt2x * kv.__bei / exf
		 << std::setw(width) << rt2x * exf * kv.__ker / _S_sqrt_pi
		 << std::setw(width) << rt2x * exf * kv.__kei / _S_sqrt_pi
		 << '\n';
	  }
	for (int i = 2001; i <= 4000; ++i)
	  {
	    auto x = _Real(0.01Q * i);
	    auto kv = std::__detail::__kelvin_asymp(nu, x);
	    auto exf = std::exp(x / _S_sqrt_2);
	    auto rt2x = _S_sqrt_2 * std::sqrt(x);
	    data << std::setw(width) << kv.__x
		 << std::setw(width) << _S_sqrt_pi * rt2x * kv.__ber / exf
		 << std::setw(width) << _S_sqrt_pi * rt2x * kv.__bei / exf
		 << std::setw(width) << rt2x * exf * kv.__ker / _S_sqrt_pi
		 << std::setw(width) << rt2x * exf * kv.__kei / _S_sqrt_pi
		 << '\n';
	  }
	data << "\n\n";
      }
  }


int
main()
{
  run_kelvin4<long double>();
  run_kelvin3<long double>();
  run_kelvin2<long double>();
  run_kelvin1<long double>();

  plot_kelvin<float>("plot/kelvin_float.txt");
  plot_kelvin<double>("plot/kelvin_double.txt");
  plot_kelvin<long double>("plot/kelvin_long_double.txt");

  plot_kelvin_order<float>("plot/kelvin_order_float.txt");
  plot_kelvin_order<double>("plot/kelvin_order_double.txt");
  plot_kelvin_order<long double>("plot/kelvin_order_long_double.txt");

  diff_kelvin2<long double>();
  diff_kelvin3<long double>();
}
