/*
g++ -std=gnu++14 -g -Wall -Wextra -I. -o test_kelvin test_kelvin.cpp
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

  /**
   * Data structure for Kelvin functions.
   */
  template<typename _Tp>
    struct _KelvinState
    {
      _Tp __nu;
      _Tp __x;
      _Tp __ber;
      _Tp __bei;
      _Tp __ker;
      _Tp __kei;
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
  template<typename _Tp>
    _Tp
    __kelvin_bex_series(_Tp __x, int __sign)
    {
      constexpr auto _S_1d2 = _Tp{1} / _Tp{2};
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __xd4 = __x / _Tp{4};
      const auto __tmp = __xd4 * __xd4;
      const auto __y = __tmp * __tmp;
      auto __term = _Tp{1};
      //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_VanWijngaardenSum<_Tp>> __bex;
      __gnu_cxx::_BasicSum<_Tp> __bex;
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
  template<typename _Tp>
    _Tp
    __kelvin_kex_series(_Tp __x, int __sign)
    {
      constexpr auto _S_1d2 = _Tp{1} / _Tp{2};
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __xd4 = __x / _Tp{4};
      const auto __tmp = __xd4 * __xd4;
      const auto __y = __tmp * __tmp;
      auto __term = _Tp{1};
      //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_VanWijngaardenSum<_Tp>> __kex;
      __gnu_cxx::_BasicSum<_Tp> __kex, _H_n;
      __kex += __term;
      _H_n += 1;
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  auto __fact = __k * (__k + __sign * _S_1d2);
	  __term *= -__y / __fact / __fact;

	  _H_n += _Tp{1} / _Tp(2 * __k);
	  auto __hterm = _Tp{1} / _Tp(2 * __k + 1);
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
  template<typename _Tp>
    _Tp
    __kelvin_ber_series(_Tp __x)
    { return __kelvin_bex_series(__x, -1); }


  /**
   * Return the Kelvin function @f$ bei(x) @f$ for real argument @c x.
   */
  template<typename _Tp>
    _Tp
    __kelvin_bei_series(_Tp __x)
    { return __x * __x * __kelvin_bex_series(__x, +1) / _Tp{4}; }


  /**
   * Return the irregular Kelvin function @f$ ker(x) @f$ for real argument @c x.
   */
  template<typename _Tp>
    _Tp
    __kelvin_ker_series(_Tp __x)
    {
      constexpr auto _S_gamma_e = __gnu_cxx::__math_constants<_Tp>::__gamma_e;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      auto __ker = __kelvin_kex_series(__x, -1);
      auto __ber = __kelvin_bex_series(__x, -1);
      auto __bei = __kelvin_bex_series(__x, +1);
      auto __xxd4 = __x * __x / _Tp{4};
      auto __ln = std::log(__x / _Tp{2}) + _S_gamma_e;
      return -__ln * __ber + _S_pi_4 * __xxd4 * __bei + __ker - _Tp{1};
    }


  /**
   * Return the irregular Kelvin function @f$ ker(x) @f$ for real argument @c x.
   */
  template<typename _Tp>
    _Tp
    __kelvin_kei_series(_Tp __x)
    {
      constexpr auto _S_gamma_e = __gnu_cxx::__math_constants<_Tp>::__gamma_e;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      auto __kei = __kelvin_kex_series(__x, +1);
      auto __ber = __kelvin_bex_series(__x, -1);
      auto __bei = __kelvin_bex_series(__x, +1);
      auto __xxd4 = __x * __x / _Tp{4};
      auto __ln = std::log(__x / _Tp{2}) + _S_gamma_e;
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
  template<typename _Tp>
    _KelvinState<_Tp>
    __kelvin_series(_Tp __x)
    {
      constexpr auto _S_gamma_e = __gnu_cxx::__math_constants<_Tp>::__gamma_e;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      constexpr auto _S_eps = _Tp{0.01L} * std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __xd2 = __x / _Tp{2};
      const auto __xxd4 = __xd2 * __xd2;
      const auto __y = __xxd4 * __xxd4;
      auto __termr = _Tp{1};
      auto __termi = _Tp{1};
      //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_VanWijngaardenSum<_Tp>>
      //  __ber, __bei, __ker, __kei;
      //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_BasicSum<_Tp>> _H_n;
      //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_BasicSum<_Tp>> __ber, __bei, __ker, __kei, _H_n;
      __gnu_cxx::_BasicSum<_Tp> __ber, __bei, __ker, __kei, _H_n;
      __ber += __termr;
      __ker += __termr;
      __bei += __termi;
      __kei += __termi;
      _H_n += _Tp{1};
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  auto __factr = _Tp{1} / _Tp(2 * __k - 1) / _Tp(2 * __k);
	  __termr *= -__y * __factr * __factr;
	  __ber += __termr;

	  _H_n += _Tp{1} / _Tp(2 * __k);
	  __ker += __termr * _H_n();

	  auto __facti = _Tp{1} / _Tp(2 * __k) / _Tp(2 * __k + 1);
	  __termi *= -__y * __facti * __facti;
	  __bei += __termi;

	  _H_n += _Tp{1} / _Tp(2 * __k + 1);
	  __kei += __termi * _H_n();

	  if (std::abs(__termr) < _S_eps * std::abs(__ber()))
	    break;
	}
      auto __ln = std::log(__x / _Tp{2}) + _S_gamma_e;
      return _KelvinState<_Tp>{_Tp{0}, __x, __ber(), __xxd4 * __bei(),
		-__ln * __ber() + _S_pi_4 * __xxd4 * __bei() + __ker() - _Tp{1},
		-__ln * __xxd4 * __bei() - _S_pi_4 * __ber() + __xxd4 * __kei()};
    }


  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename _Tp>
    _KelvinState<_Tp>
    __kelvin_asymp(_Tp __x)
    {
      constexpr auto _S_j = std::complex<_Tp>{0, 1};
      constexpr auto _S_1d2 = _Tp{1} / _Tp{2};
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      constexpr auto _S_pi_8 = _S_1d2 * _S_pi_4;
      constexpr auto _S_3pi_4 = _Tp{3} * _S_pi_4;
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Tp>::__root_pi;
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __y = _Tp{1} / (_Tp{32} * __x);
      auto __term = _Tp{1};
      const auto __xrt2 = __x / _S_sqrt_2;
      auto __barg = __xrt2 - _S_pi_8;
      auto __karg = -__xrt2 - _S_pi_8;
      //__gnu_cxx::_BasicSum<std::complex<_Tp>> __be, __ke;
      __gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_BasicSum<std::complex<_Tp>>> __be, __ke;
      __be += std::polar(__term, __barg);
      __ke += std::polar(__term, __karg);
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  __barg -= _S_pi_4;
	  __karg += _S_3pi_4;
	  auto __fact =  _Tp(2 * __k) * _Tp(2 * __k - 1) / _Tp(__k);
	  auto __next = __y * __fact * __fact / _Tp(__k);
	  if (std::abs(__next) > _Tp{1})
	    break;
	  __term *= __next;
	  __be += std::polar(__term, __barg);
	  __ke += std::polar(__term, __karg);

	  if (std::abs(__term) < _S_eps * std::abs(__be()))
	    break;
	}
      const auto __exp = std::exp(__xrt2);
      const auto __rt = std::sqrt(_Tp{2} * __x);
      const auto __kfact = _S_sqrt_pi / __rt / __exp;
      const auto __kex = __kfact * __ke();
      const auto __bfact = __exp / _S_sqrt_pi / __rt;
      const auto __bex = __bfact * __be() + _S_j * __kex / _S_pi;
      return _KelvinState<_Tp>{_Tp{0}, __x,
			       __bfact * std::real(__bex),
			       __bfact * std::imag(__bex),
			       __kfact * std::real(__kex),
			       __kfact * std::imag(__kex)};
    }


  /**
   * Compute the Kelvin functions of integer order n by series expansion.
   */
  template<typename _Tp>
    _KelvinState<_Tp>
    __kelvin_series(int __n, _Tp __x)
    {
      constexpr auto _S_j = std::complex<_Tp>{0, 1};
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      constexpr auto _S_3pi_4 = _Tp{3} * _S_pi_4;
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_maxiter = 1000;
      if (__n < 0)
	{
	  const auto _Cnp = _Tp(1 - 2 * (__n & 1));
	  auto _Kv = __kelvin_series(-__n, __x);
	  return _KelvinState<_Tp>{_Tp(__n), __x,
	  			   _Cnp * _Kv.__ber, _Cnp * _Kv.__bei,
				   _Cnp * _Kv.__ker, _Cnp * _Kv.__kei};
	}
      else
	{
	  const auto __xd2 = __x / _Tp{2};
	  const auto __y = __xd2 * __xd2;

	  __gnu_cxx::_BasicSum<std::complex<_Tp>> __be;
	  auto __bterm = _Tp{1};
	  auto __barg = _S_3pi_4 * __n;
	  __be += std::polar(__bterm, __barg);
	  for (auto __k = 1; __k < _S_maxiter; ++__k)
	    {
	      auto __bfact = (_Tp(__k) + __n) * _Tp(__k);
	      __bterm *= __y / __bfact;
	      __barg += _S_pi_2;
	      __be += std::polar(__bterm, __barg);
	      if (std::abs(__bterm) < _S_eps * std::abs(__be()))
		break;
	    }

	  __gnu_cxx::_BasicSum<std::complex<_Tp>> __ke1;
	  if (__n > 0)
	    {
	      auto __kterm1 = _Tp{1};
	      auto __karg1 = _S_3pi_4 * __n;
	      __ke1 += std::polar(__kterm1, -__karg1);
	      for (auto __k = 1; __k < __n - 1; ++__k)
		{
		  auto __fact1 = _Tp(__n - 1 - __k) * _Tp(__k);
		  __kterm1 *= __y / __fact1;
		  __karg1 += _S_pi_2;
		  __ke1 += std::polar(__kterm1, -__karg1);
		}
	      if (__n > 1)
		{
		  __kterm1 *= __y / (__n - 1);
		  __karg1 += _S_pi_2;
		  __ke1 += std::polar(__kterm1, -__karg1);
		}
	    }

	  __gnu_cxx::_BasicSum<std::complex<_Tp>> __ke2;
	  auto __hsum2 = std::__detail::__psi<_Tp>(1)
		       + std::__detail::__psi<_Tp>(1 + __n);
	  auto __kterm2 = _Tp{1};
	  auto __karg2 = _S_3pi_4 * __n;
	  __ke2 += std::polar(__hsum2 * __kterm2, __karg2);
	  for (auto __k = 1; __k < _S_maxiter; ++__k)
	    {
	      __hsum2 += _Tp{1} / _Tp(1 + __k) + _Tp{1} / _Tp(1 + __n + __k);
	      __kterm2 *= __y / _Tp(__k) /  _Tp(__n + __k);
	      __karg2 += _S_pi_2;
	      __ke2 += std::polar(__hsum2 * __kterm2, __karg2);
	      if (std::abs(__hsum2 * __kterm2) < _S_eps * std::abs(__ke2()))
		break;
	    }

	  const auto __nfact = std::__detail::__factorial<_Tp>(__n);
	  const auto __pow = std::pow(__xd2, __n);
	  auto __bex = __pow * __be() / __nfact;
	  auto __kex = -std::log(__xd2) * __bex - _S_j * _S_pi_4 * __bex
		+ __ke1() / __pow / _Tp{2}
		  / (__n > 0 ? std::__detail::__factorial<_Tp>(__n - 1) : 1)
		+ __pow * __ke2() / __nfact / _Tp{2}
		- _Tp{1}; // I needed this for ker(x) too.
	  return _KelvinState<_Tp>{_Tp(__n), __x,
				   std::real(__bex), std::imag(__bex),
				   std::real(__kex), std::imag(__kex)};
	}
    }


  /**
   * Compute the Kelvin functions of order @f$ \nu @f$ by series expansion.
   */
  template<typename _Tp>
    _KelvinState<_Tp>
    __kelvin_series(_Tp __nu, _Tp __x)
    {
      constexpr auto _S_j = std::complex<_Tp>{0, 1};
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      constexpr auto _S_3pi_4 = _Tp{3} * _S_pi_4;
      constexpr auto _S_maxiter = 1000;

      if (__gnu_cxx::__fpinteger(__nu))
	{
	  int __n = std::nearbyint(__nu);
	  return __kelvin_series(__n, __x);
	}
      else
	{
	  const auto __xd2 = __x / _Tp{2};
	  const auto __y = __xd2 * __xd2;

	  __gnu_cxx::_BasicSum<std::complex<_Tp>> __bep;
	  auto __termp = _Tp{1};
	  auto __argp = _S_3pi_4 * __nu;
	  __bep += std::polar(__termp, __argp);
	  for (auto __k = 1; __k < _S_maxiter; ++__k)
	    {
	      auto __factp = (_Tp(__k) + __nu) * _Tp(__k);
	      __termp *= __y / __factp;
	      __argp += _S_pi_2;
	      __bep += std::polar(__termp, __argp);
	    }
	  auto __bex = std::pow(__xd2, __nu) * __bep()
		     / std::tgamma(_Tp{1} + __nu);

	  __gnu_cxx::_BasicSum<std::complex<_Tp>> __bem;
	  auto __termm = _Tp{1};
	  auto __argm = -_S_3pi_4 * __nu;
	  __bem += std::polar(__termm, __argm);
	  for (auto __k = 1; __k < _S_maxiter; ++__k)
	    {
	      auto __factm = (_Tp(__k) - __nu) * _Tp(__k);
	      __termm *= __y / __factm;
	      __argm += _S_pi_2;
	      __bem += std::polar(__termm, __argm);
	    }
	  auto __bey = std::pow(__xd2, -__nu) * __bem()
		     / std::tgamma(_Tp{1} - __nu);
	  const auto __nupi = __nu * _S_pi;
	  const auto __csc = _Tp{1} / std::sin(__nupi);
	  const auto __cot = _Tp{1} / std::tan(__nupi);
	  auto __kex = _S_pi_2 * (__csc * __bey - __cot * __bex + _S_j * __bex);

	  return _KelvinState<_Tp>{__nu, __x,
				   std::real(__bex), std::imag(__bex),
				   std::real(__kex), std::imag(__kex)};
	}
    }


  /**
   * Compute the Kelvin functions by asymptotic series expansion.
   */
  template<typename _Tp>
    _KelvinState<_Tp>
    __kelvin_asymp(_Tp __nu, _Tp __x)
    {
      constexpr auto _S_j = std::complex<_Tp>{0, 1};
      constexpr auto _S_1d2 = _Tp{1} / _Tp{2};
      constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
      constexpr auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
      constexpr auto _S_pi_4 = __gnu_cxx::__math_constants<_Tp>::__pi_quarter;
      constexpr auto _S_pi_8 = _S_1d2 * _S_pi_4;
      constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
      constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Tp>::__root_pi;
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __y = _Tp{1} / (_Tp{2} * __x);
      auto __bterm = _Tp{1};
      auto __kterm = _Tp{1};
      const auto __xrt2 = __x / _S_sqrt_2;
      auto __barg = __xrt2 + __nu * _S_pi_2 - _S_pi_8;
      auto __karg = __xrt2 + __nu * _S_pi_2 + _S_pi_8;
      //__gnu_cxx::_BasicSum<std::complex<_Tp>> __be, __ke;
      __gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_BasicSum<std::complex<_Tp>>> __be, __ke;
      __be += std::polar(__bterm, __barg);
      __ke += std::polar(__kterm, __karg);
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  __barg -= _S_pi_4;
	  __karg += _S_pi_4;
	  auto __fact = (_Tp(__k) - _Tp{0.5L} - __nu)
		      * (_Tp(__k) - _Tp{0.5L} + __nu) / _Tp(__k);
	  auto __next = __y * __fact;
	  if (std::abs(__next) > _Tp{1})
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
      return _KelvinState<_Tp>{__nu, __x,
			       __bfact * std::real(__bex),
			       __bfact * std::imag(__bex),
			       __kfact * std::real(__kex),
			       -__kfact * std::imag(__kex)};
    }


/*
 *
 */
template<typename _Tp>
  void
  run_kelvin()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    __kelvin_ber_series(_Tp{});

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
	auto x = _Tp(0.1L) * i;
	auto ber = __kelvin_ber_series(x);
	auto bei = __kelvin_bei_series(x);
	auto ker = __kelvin_ker_series(x);
	auto kei = __kelvin_kei_series(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << ber
		  << std::setw(width) << bei
		  << std::setw(width) << ker
		  << std::setw(width) << kei
		  << '\n';
      }
    std::cout << std::endl;
  }


/*
 *
 */
template<typename _Tp>
  void
  run_kelvin2()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    __kelvin_series(_Tp{});

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
	auto x = _Tp(0.1L) * i;
	auto ke = __kelvin_series(x);
	std::cout << std::setw(width) << ke.__x
		  << std::setw(width) << ke.__ber
		  << std::setw(width) << ke.__bei
		  << std::setw(width) << ke.__ker
		  << std::setw(width) << ke.__kei
		  << '\n';
      }
    std::cout << std::endl;
  }


/*
 *
 */
template<typename _Tp>
  void
  diff_kelvin2()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
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
    for (int i = 200; i <= 400; ++i)
      {
	auto x = _Tp(0.1L) * i;
	auto kes = __kelvin_series(x);
	auto kea = __kelvin_asymp(x);
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


/*
 *
 */
template<typename _Tp>
  void
  run_kelvin3(_Tp nu = _Tp{0})
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
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
	auto x = _Tp(0.1L) * i;
	auto ke = __kelvin_series(nu, x);
	std::cout << std::setw(width) << ke.__x
		  << std::setw(width) << ke.__ber
		  << std::setw(width) << ke.__bei
		  << std::setw(width) << ke.__ker
		  << std::setw(width) << ke.__kei
		  << '\n';
      }
    std::cout << std::endl;
  }


/*
 *
 */
template<typename _Tp>
  void
  diff_kelvin3(_Tp nu = _Tp{0})
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
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
    for (int i = 200; i <= 400; ++i)
      {
	auto x = _Tp(0.1L) * i;
	auto kes = __kelvin_series(nu, x);
	auto kea = __kelvin_asymp(nu, x);
	std::cout << "  " << std::setw(width) << kes.__x
		  << "  " << std::setw(width) << (kea.__ber - kes.__ber) / std::abs(kes.__ber)
		  << "  " << std::setw(width) << (kea.__bei - kes.__bei) / std::abs(kes.__bei)
		  << "  " << std::setw(width) << (kea.__ker - kes.__ker) / std::abs(kes.__ker)
		  << "  " << std::setw(width) << (kea.__kei - kes.__kei) / std::abs(kes.__kei)
		  << '\n';
      }
    std::cout << std::endl;
  }


/*
 *
 */
template<typename _Tp>
  void
  run_kelvin4(int n = 0)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
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
	auto x = _Tp(0.1L) * i;
	auto ke = __kelvin_series(n, x);
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
 * 
 */
template<typename _Tp>
  void
  plot_kelvin(std::string filename)
  {
    constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
    constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Tp>::__root_pi;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Tp>::digits10);
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
	auto x = _Tp(0.01Q * i);
	auto kv = __kelvin_series(x);
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
	auto x = _Tp(0.01Q * i);
	auto kv = __kelvin_asymp(x);
	auto exf = std::exp(x / _S_sqrt_2);
	auto rt2x = _S_sqrt_2 * std::sqrt(x);
	data << std::setw(width) << kv.__x
	     << std::setw(width) << _S_sqrt_pi * _S_sqrt_pi * rt2x * rt2x * kv.__ber / exf / exf
	     << std::setw(width) << _S_sqrt_pi * _S_sqrt_pi * rt2x * rt2x * kv.__bei / exf / exf
	     << std::setw(width) << rt2x * rt2x * exf * exf * kv.__ker / _S_sqrt_pi / _S_sqrt_pi
	     << std::setw(width) << rt2x * rt2x * exf * exf * kv.__kei / _S_sqrt_pi / _S_sqrt_pi
	     << '\n';
      }
    data << "\n\n";
  }


/**
 * 
 */
template<typename _Tp>
  void
  plot_kelvin_order(std::string filename)
  {
    constexpr auto _S_sqrt_2 = __gnu_cxx::__math_constants<_Tp>::__root_2;
    constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Tp>::__root_pi;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Tp>::digits10);
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

    std::vector<_Tp> nuv{_Tp{0.0L}, _Tp{0.5L}, _Tp{1.0L}, _Tp{1.5L}, _Tp{2.0L}, _Tp{5.0L}};
    for (auto nu : nuv)
      {
	for (int i = 0; i <= 2000; ++i)
	  {
	    auto x = _Tp(0.01Q * i);
	    auto kv = __kelvin_series(nu, x);
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
	    auto x = _Tp(0.01Q * i);
	    auto kv = __kelvin_asymp(nu, x);
	    auto exf = std::exp(x / _S_sqrt_2);
	    auto rt2x = _S_sqrt_2 * std::sqrt(x);
	    data << std::setw(width) << kv.__x
		 << std::setw(width) << _S_sqrt_pi * _S_sqrt_pi * rt2x * rt2x * kv.__ber / exf / exf
		 << std::setw(width) << _S_sqrt_pi * _S_sqrt_pi * rt2x * rt2x * kv.__bei / exf / exf
		 << std::setw(width) << rt2x * rt2x * exf * exf * kv.__ker / _S_sqrt_pi / _S_sqrt_pi
		 << std::setw(width) << rt2x * rt2x * exf * exf * kv.__kei / _S_sqrt_pi / _S_sqrt_pi
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
  run_kelvin<long double>();
  plot_kelvin<double>("plot/kelvin_double.txt");
  plot_kelvin_order<double>("plot/kelvin_order_double.txt");
  diff_kelvin2<long double>();
  diff_kelvin3<long double>();
}
