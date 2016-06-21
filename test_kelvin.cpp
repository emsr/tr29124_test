/*
g++ -std=gnu++14 -g -I. -o test_kelvin test_kelvin.cpp
./test_kelvin > test_kelvin.txt
*/

#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <bits/summation.h>

bool WRITE_TERM = false;


  /*
   *
   */
  template<typename _Tp>
    _Tp
    kelvin_bex_series(_Tp __x, int __sign)
    {
      constexpr auto _S_1d2 = _Tp{1} / _Tp{2};
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __xd4 = __x / _Tp{4};
      const auto __tmp = __xd4 * __xd4;
      const auto __y = __tmp * __tmp;
      auto __fact = _Tp{1};
      auto __term = _Tp{1};
if (WRITE_TERM)
  std::cout << __term << '\n';
      //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_VanWijngaardenSum<_Tp>> __bex;
      __gnu_cxx::_BasicSum<_Tp> __bex;
      __bex += __term;
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  __fact *= __k * (__k + __sign * _S_1d2);
	  __term *= -__y / __fact / __fact;
if (WRITE_TERM)
  std::cout << __term << '\n';
	  __bex += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__bex()))
	    break;
	}
      return __bex();
    }


  /*
   *
   */
  template<typename _Tp>
    _Tp
    kelvin_ber_series(_Tp __x)
    { return kelvin_bex_series(__x, -1); }


  /*
   *
   */
  template<typename _Tp>
    _Tp
    kelvin_bei_series(_Tp __x)
    { return __x * __x * kelvin_bex_series(__x, +1) / _Tp{4}; }


  /*
   *
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    kelvin_series(_Tp __x)
    {
      constexpr auto _S_gamma_e = __gnu_cxx::__math_constants<_Tp>::__gamma_e;
      constexpr auto _S_1d2 = _Tp{1} / _Tp{2};
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __xd2 = __x / _Tp{2};
      const auto __xxd4 = __xd2 * __xd2;
      const auto __y = __xxd4 * __xxd4;
      auto __factr = _Tp{1};
      auto __termr = _Tp{1};
      auto __facti = _Tp{1};
      auto __termi = _Tp{1};
if (WRITE_TERM)
  std::cout << __termr << '\t' << __termi << '\n';
      //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_VanWijngaardenSum<_Tp>>
      //  __ber, __bei, __ker, __kei;
      //__gnu_cxx::_WenigerDeltaSum<__gnu_cxx::_BasicSum<_Tp>> _H_n;
      __gnu_cxx::_BasicSum<_Tp> __ber, __bei, __ker, __kei;
      __ber += __termr;
      __bei += __termi;
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  __factr *= _Tp{2} * __k * (_Tp{2} * __k - 1);
	  __termr *= -__y / __factr / __factr;
	  __ber += __termr;

	  __facti = __factr * (_Tp{2} * __k + 1) / (_Tp{2} * __k - 1);
	  __termi *= -__y / __facti / __facti;
	  __bei += __termi;

if (WRITE_TERM)
  std::cout << __termr << '\t' << __termi << '\n';
	  if (std::abs(__termr) < _S_eps * std::abs(__ber()))
	    break;
	}
      return std::make_pair(__ber(), __xxd4 * __bei());
    }


  /*
   *
   */
  template<typename _Tp>
    std::pair<_Tp, _Tp>
    kelvin_series2(_Tp __x)
    {
      constexpr auto _S_gamma_e = __gnu_cxx::__math_constants<_Tp>::__gamma_e;
      constexpr auto _S_1d2 = _Tp{1} / _Tp{2};
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr auto _S_maxiter = 1000;
      const auto __xd2 = __x / _Tp{2};
      const auto __xxd4 = __xd2 * __xd2;
      const auto __y = __xxd4 * __xxd4;
      auto __factr = _Tp{1};
      auto __termr = _Tp{1};
      auto __ber = __termr;
      auto __facti = _Tp{1};
      auto __termi = _Tp{1};
      auto __bei = __termi;
if (WRITE_TERM)
  std::cout << __termr << '\t' << __termi << '\n';
      __ber += __termr;
      __bei += __termi;
      for (auto __k = 1; __k < _S_maxiter; ++__k)
	{
	  __factr *= _Tp{2} * __k * (_Tp{2} * __k - 1);
	  __termr *= -__y / __factr / __factr;
	  __ber += __termr;

	  __facti = __factr * (_Tp{2} * __k + 1) / (_Tp{2} * __k - 1);
	  __termi *= -__y / __facti / __facti;
	  __bei += __termi;

if (WRITE_TERM)
  std::cout << __termr << '\t' << __termi << '\n';
	  if (std::abs(__termr) < _S_eps * std::abs(__ber))
	    break;
	}
      return std::make_pair(__ber, __xxd4 * __bei);
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

    auto ber0 = kelvin_ber_series(_Tp{});

WRITE_TERM=true;
    auto ber1 = kelvin_ber_series(_Tp{1});
WRITE_TERM=false;

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
    for (int i = -100; i <= 200; ++i)
      {
	auto x = _Tp(0.1) * i;
	auto ber = kelvin_ber_series(x);
	auto bei = kelvin_bei_series(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << ber
		  << std::setw(width) << bei
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

    auto ber0 = kelvin_series2(_Tp{});

WRITE_TERM=true;
    auto ber1 = kelvin_series2(_Tp{1});
WRITE_TERM=false;

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
    for (int i = -100; i <= 200; ++i)
      {
	auto x = _Tp(0.1) * i;
	auto be = kelvin_series2(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << be.first
		  << std::setw(width) << be.second
		  << '\n';
      }
    std::cout << std::endl;
  }


int
main()
{
  run_kelvin2<double>();
}
