
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

#include <emsr/continued_fractions.h>

template<typename Tp>
  void
  test_lentz_continued_fraction(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto S_pi_2 = emsr::pi_v<Tp> / Tp{2};
    using Cmplx = std::complex<Tp>;

    auto a_trigint
      = [](std::size_t i, Tp)
	-> Cmplx
	{
	  if (i == 1)
	    return Tp(1);
	  else
	    return -Tp(i - 1) * Tp(i - 1);
	};
    using AFun = decltype(a_trigint);

    auto b_trigint
      = [](std::size_t i, Tp x)
	{
	  if (i == 0)
	    return Cmplx{0};
	  else
	    return Cmplx{Tp(2 * i - 1), x};
	};
    using BFun = decltype(b_trigint);

    auto w_trigint
      = [](std::size_t, Tp)
	{ return Cmplx{}; };
    using TailFun = decltype(w_trigint);

    emsr::LentzContinuedFraction<Tp, AFun, BFun, TailFun>
      SiCi(a_trigint, b_trigint, w_trigint);
    auto t = 1.2;
    auto y = SiCi(t);
    y *= std::polar(Tp{1}, -t);
    y += Cmplx{0, S_pi_2};
    std::cout << '\n';
    std::cout << "SiCi = " << y << '\n';
    std::cout << '\n';
    std::cout << "Si = " << std::setw(width) << std::imag(y) << '\n';
    std::cout << "Ci = " << std::setw(width) << -std::real(y) << '\n';
    std::cout << '\n';
    std::cout << "Si = " << std::setw(width) << __gnu_cxx::sinint(1.2) << '\n';
    std::cout << "Ci = " << std::setw(width) << __gnu_cxx::cosint(1.2) << '\n';
  }

template<typename Tp>
  void
  test_lentz_hypint(Tp proto = Tp{})
  {
    using Val = Tp;
    using Real = emsr::num_traits_t<Val>;
    using Cmplx = std::complex<Real>;

    std::cout.precision(emsr::digits10(std::real(proto)));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto S_pi_2 = emsr::pi_v<Tp> / Tp{2};
    const auto S_i = Cmplx{0, 1};

    auto a_hypint
      = [](std::size_t i, Tp)
	-> Cmplx
	{
	  if (i == 1)
	    return Tp(1);
	  else
	    return -Tp(i - 1) * Tp(i - 1);
	};
    using AFun = decltype(a_hypint);

    auto b_hypint
      = [S_i](std::size_t i, Tp x)
	{
	  if (i == 0)
	    return Cmplx{0};
	  else
	    return Cmplx(Tp(2 * i - 1), x);
	};
    using BFun = decltype(b_hypint);

    auto w_hypint
      = [](std::size_t, Tp)
	{ return Cmplx{}; };
    using TailFun = decltype(w_hypint);

    emsr::LentzContinuedFraction<Tp, AFun, BFun, TailFun>
      ShiChi(a_hypint, b_hypint, w_hypint);
    auto t = 1.2;
    auto y = ShiChi(t);
    y *= std::polar(Tp{1}, -t);
    y += Cmplx{0, S_pi_2};
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
