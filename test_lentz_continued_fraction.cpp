/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -I. -o test_lentz_continued_fraction test_lentz_continued_fraction.cpp -lquadmath
./test_lentz_continued_fraction > test_lentz_continued_fraction.txt

$HOME/bin/bin/g++ -std=gnu++17 -I. -o test_lentz_continued_fraction test_lentz_continued_fraction.cpp -lquadmath
./test_lentz_continued_fraction > test_lentz_continued_fraction.txt
*/

#include <ext/cmath>
#include <complex>
#include <iostream>

#include "LentzContinuedFraction.tcc"

template<typename _Tp>
  void
  test_lentz_continued_fraction(_Tp proto = _Tp{})
  {
    const auto _S_pi_2 = __gnu_cxx::__const_pi_half(proto);
    using _Cmplx = std::complex<_Tp>;

    auto a_trigint
      = [](size_t i, _Tp)
        -> _Cmplx
        { return -_Tp(i) * _Tp(i); };
    using _AFun = decltype(a_trigint);

    auto b_trigint
      = [](size_t i, _Tp __x)
	{
	  if (i == 0)
	    return _Cmplx{0};
	  else if (i == 1)
	    return _Cmplx{1, __x};
	  else
	    return _Cmplx{1, __x} + _Tp(i * 2);
	};
    using _BFun = decltype(b_trigint);

    auto w_trigint
      = [](size_t, _Tp)
	{ return _Cmplx{}; };
    using _TailFun = decltype(w_trigint);

    _LentzContinuedFraction<_Tp, _AFun, _BFun, _TailFun>
      SiCi(a_trigint, b_trigint, w_trigint);
    auto t = 1.2;
    auto y = _Cmplx{0, _S_pi_2} + SiCi(t);
    y *= std::polar(_Tp{1}, -t);
    std::cout << "SiCi = " << y << '\n';
    std::cout << "Si = " << std::imag(y) << '\n';
    std::cout << "Ci = " << -std::real(y) << '\n';
    std::cout << "Si = " << __gnu_cxx::sinint(1.2) << '\n';
    std::cout << "Ci = " << __gnu_cxx::cosint(1.2) << '\n';
  }

int
main()
{
  test_lentz_continued_fraction(1.0);
}
