
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

#include <ext/continued_fractions.h>

template<typename _Tp>
  void
  test_lentz_continued_fraction(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto _S_pi_2 = emsr::pi_v<_Tp> / _Tp{2};
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

    _LentzContinuedFraction<_Tp, _AFun, _BFun, _TailFun>
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
  test_lentz_hypint(_Tp proto = _Tp{})
  {
    using _Val = _Tp;
    using _Real = emsr::num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    std::cout.precision(__gnu_cxx::__digits10(std::real(proto)));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto _S_pi_2 = emsr::pi_v<_Tp> / _Tp{2};
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

    _LentzContinuedFraction<_Tp, _AFun, _BFun, _TailFun>
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
  test_lentz_continued_fraction(1.0);
  test_lentz_hypint(1.0);
}
