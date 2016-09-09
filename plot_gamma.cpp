/*
$HOME/bin_tr29124/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -I. -o plot_gamma plot_gamma.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./plot_gamma > plot_gamma.txt

$HOME/bin/bin/g++ -std=gnu++14 -DNO_LOGBQ -I. -o plot_gamma plot_gamma.cpp -lquadmath
./plot_gamma > plot_gamma.txt
*/

#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <bits/specfun.h>


template<typename _Tp>
  void
  plot_spouge(std::string filename)
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    using GammaT = decltype(std::__detail::__log_gamma_spouge(_Cmplx{}));
    std::vector<std::vector<GammaT>> gammav;

    for (int i = -200; i <= +50; ++i)
      {
	gammav.push_back(std::vector<GammaT>{});
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Cmplx(0.10Q * i, 0.10Q * j);
	    gammav.back().push_back(std::__detail::__log_gamma_spouge(t));
	  }
      }

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto z = _Cmplx(0.10Q * i, 0.10Q * j);
	    auto gamma = gammav[i][j];
	    data << std::setw(width) << std::real(z)
		 << std::setw(width) << std::imag(z)
		 << std::setw(width) << std::real(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto z = _Cmplx(0.10Q * i, 0.10Q * j);
	    auto gamma = gammav[i][j];
	    data << std::setw(width) << std::real(z)
		 << std::setw(width) << std::imag(z)
		 << std::setw(width) << std::imag(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto z = _Cmplx(0.10Q * i, 0.10Q * j);
	    auto gamma = gammav[i][j];
	    data << std::setw(width) << std::real(z)
		 << std::setw(width) << std::imag(z)
		 << std::setw(width) << std::abs(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto z = _Cmplx(0.10Q * i, 0.10Q * j);
	    auto gamma = gammav[i][j];
	    data << std::setw(width) << std::real(z)
		 << std::setw(width) << std::imag(z)
		 << std::setw(width) << std::arg(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';
  }

template<typename _Tp>
  void
  plot_lanczos(std::string filename)
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    using GammaT = decltype(std::__detail::__log_gamma_lanczos(_Cmplx{}));
    std::vector<std::vector<GammaT>> gammav;

    for (int i = -200; i <= +50; ++i)
      {
	gammav.push_back(std::vector<GammaT>{});
	for (int j = -50; j <= +50; ++j)
	  {
	    auto t = _Cmplx(0.10Q * i, 0.10Q * j);
	    gammav.back().push_back(std::__detail::__log_gamma_lanczos(t));
	  }
      }

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto z = _Cmplx(0.10Q * i, 0.10Q * j);
	    auto gamma = gammav[i][j];
	    data << std::setw(width) << std::real(z)
		 << std::setw(width) << std::imag(z)
		 << std::setw(width) << std::real(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto z = _Cmplx(0.10Q * i, 0.10Q * j);
	    auto gamma = gammav[i][j];
	    data << std::setw(width) << std::real(z)
		 << std::setw(width) << std::imag(z)
		 << std::setw(width) << std::imag(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto z = _Cmplx(0.10Q * i, 0.10Q * j);
	    auto gamma = gammav[i][j];
	    data << std::setw(width) << std::real(z)
		 << std::setw(width) << std::imag(z)
		 << std::setw(width) << std::abs(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = -200; i <= +50; ++i)
      {
	for (int j = -50; j <= +50; ++j)
	  {
	    auto z = _Cmplx(0.10Q * i, 0.10Q * j);
	    auto gamma = gammav[i][j];
	    data << std::setw(width) << std::real(z)
		 << std::setw(width) << std::imag(z)
		 << std::setw(width) << std::arg(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';
  }

int
main()
{
  std::cout << "\n\nLanczos Algorithm\n\n";
  std::cout << "\nlanczos<float>\n";
  plot_lanczos<float>("plot/gamma_lanczos_float.txt");
  std::cout << "\nlanczos<double>\n";
  plot_lanczos<double>("plot/gamma_lanczos_double.txt");
  std::cout << "\nlanczos<long double>\n";
  plot_lanczos<long double>("plot/gamma_lanczos_long_double.txt");
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\nlanczos<__float128>\n";
  plot_lanczos<__float128>("plot/gamma_lanczos__float128.txt");
#endif

  std::cout << "\n\nSpouge Algorithm\n\n";
  std::cout << "\nspouge<float>\n";
  plot_spouge<float>("plot/gamma_spouge_float.txt");
  std::cout << "\nspouge<double>\n";
  plot_spouge<double>("plot/gamma_spouge_double.txt");
  std::cout << "\nspouge<long double>\n";
  plot_spouge<long double>("plot/gamma_spouge_long_double.txt");
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\nspouge<__float128>\n";
  plot_spouge<__float128>("plot/gamma_spouge__float128.txt");
#endif
}
