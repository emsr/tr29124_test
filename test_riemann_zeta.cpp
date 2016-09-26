/*
$HOME/bin_tr29124/bin/g++ -std=c++17 -g -I. -o test_riemann_zeta test_riemann_zeta.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_riemann_zeta > test_riemann_zeta.txt

$HOME/bin/bin/g++ -std=gnu++14 -DNO_LOGBQ -I. -o test_riemann_zeta test_riemann_zeta.cpp -lquadmath
./test_riemann_zeta > test_riemann_zeta.txt
*/

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

// I'm not sure why I need this here and not other places...
template<>
  constexpr std::array<float, 7>
  std::__detail::_GammaLanczos<float>::_S_cheby;
template<>
  constexpr std::array<double, 10>
  std::__detail::_GammaLanczos<double>::_S_cheby;
template<>
  constexpr std::array<long double, 11>
  std::__detail::_GammaLanczos<long double>::_S_cheby;

template<typename _Tp>
  void
  plot_riemann_zeta(std::string filename)
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    constexpr auto deg = __gnu_cxx::__math_constants<_Real>::__deg;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    using zetaT = decltype(std::__detail::__riemann_zeta(_Cmplx{}));
    std::vector<std::vector<zetaT>> sv;
    std::vector<std::vector<zetaT>> zetav;

    int i_min = -200;
    int j_min = -50;

    for (int i = i_min; i <= +50; ++i)
      {
        sv.push_back(std::vector<zetaT>{});
	zetav.push_back(std::vector<zetaT>{});
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = _Cmplx(0.10L * i, 0.10L * j);
	    sv.back().push_back(s);
	    zetav.back().push_back(std::__detail::__riemann_zeta(s));
	  }
      }

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = sv[i - i_min][j - j_min];
	    auto zeta = zetav[i - i_min][j - j_min];
	    data << std::setw(width) << std::real(s)
		 << std::setw(width) << std::imag(s)
		 << std::setw(width) << std::real(zeta)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = sv[i - i_min][j - j_min];
	    auto zeta = zetav[i - i_min][j - j_min];
	    data << std::setw(width) << std::real(s)
		 << std::setw(width) << std::imag(s)
		 << std::setw(width) << std::imag(zeta)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = sv[i - i_min][j - j_min];
	    auto zeta = zetav[i - i_min][j - j_min];
	    data << std::setw(width) << std::real(s)
		 << std::setw(width) << std::imag(s)
		 << std::setw(width) << std::abs(zeta)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = sv[i - i_min][j - j_min];
	    auto zeta = zetav[i - i_min][j - j_min];
	    data << std::setw(width) << std::real(s)
		 << std::setw(width) << std::imag(s)
		 << std::setw(width) << deg * std::arg(zeta) 
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';
  }

template<typename _Tp>
  void
  test_riemann_zeta_real()
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    constexpr auto deg = __gnu_cxx::__math_constants<_Real>::__deg;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    using zetaT = decltype(std::__detail::__riemann_zeta(_Cmplx{}));
    std::vector<std::vector<zetaT>> sv;
    std::vector<std::vector<zetaT>> zetav;

    int i_min = -250;

    std::cout << '\n'
	      << std::setw(width) << "s"
	      << std::setw(4 + 2 * width) << "zetac = zeta (cmplx)"
	      << std::setw(width) << "zeta (real)"
	      << std::setw(width) << "|zetac - zeta|"
	      << '\n';
    auto ac = _Cmplx(1.0L);
    auto a = _Real(1.0L);
    for (int i = i_min; i <= +250; ++i)
      {
        auto sc = _Cmplx(0.10L * i, 0.0L);
	auto zetac = std::__detail::__riemann_zeta(sc);
        auto s = _Real(0.10L * i);
	auto zeta = std::__detail::__riemann_zeta(s);
	std::cout << std::setw(width) << s
		  << std::setw(4 + 2 * width) << zetac
		  << std::setw(width) << zeta
		  << std::setw(width) << std::abs(zetac - zeta)
		  << '\n';
      }
  }

int
main()
{
  using namespace std::literals::complex_literals;

  auto zetam = std::__detail::__riemann_zeta(0.01l - 1.0il);
  std::cout << "zeta(" << 0.01l - 1.0il << ") = " << zetam << '\n';
  auto zetap = std::__detail::__riemann_zeta(0.01l + 1.0il);
  std::cout << "zeta(" << 0.01l + 1.0il << ") = " << zetap << '\n';

  test_riemann_zeta_real<long double>();

  std::cout << "\n\nRiemann zeta\n\n";

  std::cout << "\nriemann_zeta<float>\n";
  plot_riemann_zeta<float>("plot/riemann_zeta_float.txt");

  std::cout << "\nriemann_zeta<double>\n";
  plot_riemann_zeta<double>("plot/riemann_zeta_double.txt");

  std::cout << "\nriemann_zeta<long double>\n";
  plot_riemann_zeta<long double>("plot/riemann_zeta_long_double.txt");

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //std::cout << "\nriemann_zeta<__float128>\n";
  //plot_riemann_zeta<__float128>("plot/riemann_zeta__float128.txt");
#endif
}
