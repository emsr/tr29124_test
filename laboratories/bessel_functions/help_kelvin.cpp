/*
$HOME/bin/bin/g++ -std=c++14 -I../../include -o help_kelvin help_kelvin.cpp
./help_kelvin > help_kelvin.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <ext/math_const.h>

template<typename _Tp>
  std::pair<_Tp, _Tp>
  kelvin_ber(_Tp __x)
  {
    constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    constexpr int _S_max_iter = 100;

    if (__x == _Tp{0})
      return std::make_pair(_Tp{1}, _Tp{0});

    const auto __y = __x / _Tp{2};
    const auto __y2 = __y * __y;
    const auto __y4 = __y2 * __y2;

    auto __term = _Tp{1};
    auto __ber = __term;
    auto __berd = _Tp{0};
    for (int __k = 1; __k < _S_max_iter; ++__k)
      {
	const auto __fact = _Tp{1} / _Tp(2 * __k - 1) / _Tp(2 * __k);
	__term *= -__y4 * __fact * __fact;
	__ber += __term;
	__berd += _Tp(2 * __k) * __term / __y;
	if (std::abs(__term) < _S_eps * std::abs(__ber))
	  break;
      }

    return std::make_pair(__ber, __berd);
  }

template<typename _Tp>
  std::pair<_Tp, _Tp>
  kelvin_bei(_Tp __x)
  {
    constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    constexpr int _S_max_iter = 100;

    if (__x == _Tp{0})
      return std::make_pair(_Tp{0}, _Tp{0});

    const auto __y = __x / _Tp{2};
    const auto __y2 = __y * __y;
    const auto __y4 = __y2 * __y2;

    auto __term = __y2;
    auto __bei = __term;
    auto __beid = __y;
    for (int __k = 1; __k < _S_max_iter; ++__k)
      {
	const auto __fact = _Tp{1} / _Tp(2 * __k + 1) / _Tp(2 * __k);
	__term *= -__y4 * __fact * __fact;
	__bei += __term;
	__beid += _Tp(2 * __k + 1) * __term / __y;
	if (std::abs(__term) < _S_eps * std::abs(__bei))
	  break;
      }

    return std::make_pair(__bei, __beid);
  }

template<typename _Tp>
  std::pair<_Tp, _Tp>
  kelvin_ker(_Tp __x)
  {
    constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    constexpr int _S_max_iter = 100;
    constexpr auto _S_gamma_e = __gnu_cxx::__const_gamma_e<_Tp>();
    constexpr auto _S_pi_4 = __gnu_cxx::__const_pi_quarter<_Tp>();
    constexpr auto _S_inf = std::numeric_limits<_Tp>::infinity();

    if (__x == _Tp{0})
      return std::make_pair(_S_inf, -_S_inf);

    const auto __y = __x / _Tp{2};
    const auto __y2 = __y * __y;
    const auto __y4 = __y2 * __y2;

    // phi(1) = -gamma_E <=> H_0 = 0.
    auto _H = _Tp{0};
    auto __term = _Tp{1};
    auto __ker = _Tp{0}; 
    auto __kerd = _Tp{0};
    for (int __k = 1; __k < _S_max_iter; ++__k)
      {
	_H += _Tp{1} / _Tp(2 * __k - 1) + _Tp{1} / _Tp(2 * __k);
	const auto __fact = _Tp{1} / _Tp(2 * __k - 1) / _Tp(2 * __k);
	__term *= -__y4 * __fact * __fact;
	__ker += __term * _H;
	__kerd += _Tp(2 * __k) * __term * _H / __y;
	if (std::abs(__term * _H) < _S_eps * std::abs(__ker))
	  break;
      }

    const auto __ber = kelvin_ber(__x);
    const auto __bei = kelvin_bei(__x);
    __ker = __ker
	  - (std::log(__y) + _S_gamma_e) * __ber.first
	  + _S_pi_4 * __bei.first;
    __kerd = __kerd
	   - __ber.first / __x - (std::log(__y) + _S_gamma_e) * __ber.second
	   + _S_pi_4 * __bei.second;

    return std::make_pair(__ker, __kerd);
  }

template<typename _Tp>
  std::pair<_Tp, _Tp>
  kelvin_kei(_Tp __x)
  {
    constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    constexpr int _S_max_iter = 100;
    constexpr auto _S_gamma_e = __gnu_cxx::__const_gamma_e<_Tp>();
    constexpr auto _S_pi_4 = __gnu_cxx::__const_pi_quarter<_Tp>();

    if (__x == _Tp{0})
      return std::make_pair(-_S_pi_4, _Tp{0});

    const auto __y = __x / _Tp{2};
    const auto __y2 = __y * __y;
    const auto __y4 = __y2 * __y2;

    auto _H = _Tp{1};
    auto __term = __y2;
    auto __kei = __term;
    auto __keid = __y;
    for (int __k = 1; __k < _S_max_iter; ++__k)
      {
	_H += _Tp{1} / _Tp(2 * __k) + _Tp{1} / _Tp(2 * __k + 1);
	const auto __fact = _Tp{1} / _Tp(2 * __k) / _Tp(2 * __k + 1);
	__term *= -__y4 * __fact * __fact;
	__kei += __term * _H;
	__keid += _Tp(2 * __k + 1) * __term * _H / __y;
	if (std::abs(__term * _H) < _S_eps * std::abs(__kei))
	  break;
      }

    const auto __ber = kelvin_ber(__x);
    const auto __bei = kelvin_bei(__x);
    __kei = __kei
	  - (std::log(__y) + _S_gamma_e) * __bei.first
	  - _S_pi_4 * __ber.first;
    __keid = __keid
	   - __bei.first / __x - (std::log(__y) + _S_gamma_e) * __bei.second
	   - _S_pi_4 * __ber.second;

    return std::make_pair(__kei, __keid);
  }

template<typename _Tp>
  void
  help_kelvin()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    std::cout << "\n\nPrint Kelvin functions computed by series expansions\n";
    std::cout << std::setw(w) << "x"
	      << std::setw(w) << "ber(x)"
	      << std::setw(w) << "bei(x)"
	      << std::setw(w) << "ker(x)"
	      << std::setw(w) << "kei(x)"
	      << std::setw(w) << "ber'(x)"
	      << std::setw(w) << "bei'(x)"
	      << std::setw(w) << "ker'(x)"
	      << std::setw(w) << "kei'(x)"
	      << '\n';
    const auto del = _Tp{1}/_Tp{10};
    for (int i = 0; i <= 200; ++i)
      {
	const auto x = del * i;
	const auto ber = kelvin_ber(x);
	const auto bei = kelvin_bei(x);
	const auto ker = kelvin_ker(x);
	const auto kei = kelvin_kei(x);
	std::cout << std::setw(w) << x
		  << std::setw(w) << ber.first
		  << std::setw(w) << bei.first
		  << std::setw(w) << ker.first
		  << std::setw(w) << kei.first
		  << std::setw(w) << ber.second
		  << std::setw(w) << bei.second
		  << std::setw(w) << ker.second
		  << std::setw(w) << kei.second
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }

int
main()
{
  help_kelvin<double>();
}
