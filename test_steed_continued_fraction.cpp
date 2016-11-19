/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_steed_continued_fraction test_steed_continued_fraction.cpp -lquadmath
./test_steed_continued_fraction > test_steed_continued_fraction.txt

$HOME/bin/bin/g++ -std=gnu++17 -I. -o test_steed_continued_fraction test_steed_continued_fraction.cpp -lquadmath
./test_steed_continued_fraction > test_steed_continued_fraction.txt
*/

#include <ext/cmath>
#include <complex>
#include <iostream>
#include <iomanip>

#include "SteedContinuedFraction.tcc"

template<typename _Tp>
  void
  test_steed_continued_fraction(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto _S_pi_2 = __gnu_cxx::__const_pi_half(proto);
    using _Cmplx = std::complex<_Tp>;

    auto a_trigint
      = [](std::size_t i, _Tp)
	-> _Cmplx
	{
	  if (i == 1)
	    return _Tp(1);
	  else
	    return -_Tp(i - 1) * _Tp(i - 1);
	};
    using _AFun = decltype(a_trigint);

    auto b_trigint
      = [](std::size_t i, _Tp __x)
	{
	  if (i == 0)
	    return _Cmplx{0};
	  else
	    return _Cmplx{_Tp(2 * i - 1), __x};
	};
    using _BFun = decltype(b_trigint);

    auto w_trigint
      = [](std::size_t, _Tp)
	{ return _Cmplx{}; };
    using _TailFun = decltype(w_trigint);

    _SteedContinuedFraction<_Tp, _AFun, _BFun, _TailFun>
      SiCi(a_trigint, b_trigint, w_trigint);
    auto t = 1.2;
    auto y = SiCi(t);
    y *= std::polar(_Tp{1}, -t);
    y += _Cmplx{0, _S_pi_2};
    std::cout << '\n';
    std::cout << "SiCi = " << y << '\n';
    std::cout << '\n';
    std::cout << "Si = " << std::setw(width) << std::imag(y) << '\n';
    std::cout << "Ci = " << std::setw(width) << -std::real(y) << '\n';
    std::cout << '\n';
    std::cout << "Si = " << std::setw(width) << __gnu_cxx::sinint(1.2) << '\n';
    std::cout << "Ci = " << std::setw(width) << __gnu_cxx::cosint(1.2) << '\n';
  }

template<typename _Tp>
  void
  test_steed_hypint(_Tp proto = _Tp{})
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    std::cout.precision(__gnu_cxx::__digits10(std::real(proto)));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto _S_pi_2 = __gnu_cxx::__const_pi_half(proto);
    const auto _S_i = _Cmplx{0, 1};

    auto a_hypint
      = [](std::size_t i, _Tp)
	-> _Cmplx
	{
	  if (i == 1)
	    return _Tp(1);
	  else
	    return -_Tp(i - 1) * _Tp(i - 1);
	};
    using _AFun = decltype(a_hypint);

    auto b_hypint
      = [_S_i](std::size_t i, _Tp __x)
	{
	  if (i == 0)
	    return _Cmplx{0};
	  else
	    return _Cmplx(_Tp(2 * i - 1), __x);
	};
    using _BFun = decltype(b_hypint);

    auto w_hypint
      = [](std::size_t, _Tp)
	{ return _Cmplx{}; };
    using _TailFun = decltype(w_hypint);

    _SteedContinuedFraction<_Tp, _AFun, _BFun, _TailFun>
      ShiChi(a_hypint, b_hypint, w_hypint);
    auto t = 1.2;
    auto y = ShiChi(t);
    y *= std::polar(_Tp{1}, -t);
    y += _Cmplx{0, _S_pi_2};
    std::cout << '\n';
    std::cout << "ShiChi = " << y << '\n';
    std::cout << '\n';
    std::cout << "Shi = " << std::setw(width) << std::imag(y) << '\n';
    std::cout << "Chi = " << std::setw(width) << -std::real(y) << '\n';
    std::cout << '\n';
    std::cout << "Shi = " << std::setw(width) << __gnu_cxx::sinhint(1.2) << '\n';
    std::cout << "Chi = " << std::setw(width) << __gnu_cxx::coshint(1.2) << '\n';
  }


int
main()
{
  test_steed_continued_fraction(1.0);
  test_steed_hypint(1.0);
}
