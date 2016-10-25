/*
$HOME/bin_tr29124/bin/g++ -std=c++17 -g -I. -o test_riemann_zeta test_riemann_zeta.cpp -lquadmath
./test_riemann_zeta > test_riemann_zeta.txt

$HOME/bin/bin/g++ -std=gnu++17 -DNO_LOGBQ -I. -o test_riemann_zeta test_riemann_zeta.cpp -lquadmath
./test_riemann_zeta > test_riemann_zeta.txt
*/

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <bits/float128_io.h>

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

template<typename _Tp>
  void
  test_nontrivial_zeros()
  {
    using namespace std::literals::complex_literals;

    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    std::cout.precision(std::numeric_limits<_Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    _Cmplx
    __zeros[10]
    {
      0.5l + 14.134725141734693790457251983562470270784257115699il,
      0.5l + 21.022039638771554992628479593896902777334340524903il,
      0.5l + 25.010857580145688763213790992562821818659549672558il,
      0.5l + 30.424876125859513210311897530584091320181560023715il,
      0.5l + 32.935061587739189690662368964074903488812715603517il,
      0.5l + 37.586178158825671257217763480705332821405597350831il,
      0.5l + 40.918719012147495187398126914633254395726165962777il,
      0.5l + 43.327073280914999519496122165406805782645668371837il,
      0.5l + 48.005150881167159727942472749427516041686844001144il,
      0.5l + 49.773832477672302181916784678563724057723178299677il
    };

    std::cout << '\n'
	      << std::setw(4 + 2 * width) << "s"
	      << std::setw(4 + 2 * width) << "zeta(s)"
	      << std::setw(width) << "|zeta(s)|"
	      << '\n';
    for (auto s : __zeros)
      {
	auto zeta = std::__detail::__riemann_zeta(s);
	std::cout << std::setw(4 + 2 * width) << s
		  << std::setw(4 + 2 * width) << zeta
		  << std::setw(width) << std::abs(zeta)
		  << '\n';
      }
  }

int
main()
{
  using namespace std::literals::complex_literals;

  // These barf on Cygwin because the literals get turned to __complex__ long double!
  //auto zetam = std::__detail::__riemann_zeta(0.01l - 1.0il);
  //std::cout << "zeta(" << 0.01l - 1.0il << ") = " << zetam << '\n';
  //auto zetap = std::__detail::__riemann_zeta(0.01l + 1.0il);
  //std::cout << "zeta(" << 0.01l + 1.0il << ") = " << zetap << '\n';

  test_nontrivial_zeros<long double>();

  test_riemann_zeta_real<long double>();

  std::cout << "\n\nRiemann zeta\n\n";

  std::cout << "\nriemann_zeta<float>\n";
  plot_riemann_zeta<float>("plot/riemann_zeta_float.txt");

  std::cout << "\nriemann_zeta<double>\n";
  plot_riemann_zeta<double>("plot/riemann_zeta_double.txt");

  std::cout << "\nriemann_zeta<long double>\n";
  plot_riemann_zeta<long double>("plot/riemann_zeta_long_double.txt");

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\nriemann_zeta<__float128>\n";
  plot_riemann_zeta<__float128>("plot/riemann_zeta__float128.txt");
#endif
}
