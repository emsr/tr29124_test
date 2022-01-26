/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

template<typename _Tp>
  void
  plot_spouge(std::string filename)
  {
    using _Val = _Tp;
    using _Real = emsr::num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    constexpr auto deg = emsr::deg_v<_Real>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    using emsr::detail::spouge_log_gamma1p;
    using GammaT = decltype(spouge_log_gamma1p(_Cmplx{}));
    std::vector<std::vector<GammaT>> zv;
    std::vector<std::vector<GammaT>> gammav;

    int i_min = -200;
    int j_min = -50;

    for (int i = i_min; i <= +50; ++i)
      {
        zv.push_back(std::vector<GammaT>{});
	gammav.push_back(std::vector<GammaT>{});
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto t = _Cmplx(0.10L * i, 0.10L * j);
	    zv.back().push_back(t);
	    gammav.back().push_back(spouge_log_gamma1p(t - GammaT{1}));
	  }
      }

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto z = zv[i - i_min][j - j_min];
	    auto gamma = gammav[i - i_min][j - j_min];
	    data << std::setw(w) << std::real(z)
		 << std::setw(w) << std::imag(z)
		 << std::setw(w) << std::real(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto z = zv[i - i_min][j - j_min];
	    auto gamma = gammav[i - i_min][j - j_min];
	    data << std::setw(w) << std::real(z)
		 << std::setw(w) << std::imag(z)
		 << std::setw(w) << std::imag(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto z = zv[i - i_min][j - j_min];
	    auto gamma = gammav[i - i_min][j - j_min];
	    data << std::setw(w) << std::real(z)
		 << std::setw(w) << std::imag(z)
		 << std::setw(w) << std::abs(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto z = zv[i - i_min][j - j_min];
	    auto gamma = gammav[i - i_min][j - j_min];
	    data << std::setw(w) << std::real(z)
		 << std::setw(w) << std::imag(z)
		 << std::setw(w) << deg * std::arg(gamma) 
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
    using _Real = emsr::num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    constexpr auto deg = emsr::deg_v<_Real>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    using emsr::detail::lanczos_log_gamma1p;
    using GammaT = decltype(lanczos_log_gamma1p(_Cmplx{}));
    std::vector<std::vector<GammaT>> zv;
    std::vector<std::vector<GammaT>> gammav;

    int i_min = -200;
    int j_min = -50;

    for (int i = i_min; i <= +50; ++i)
      {
        zv.push_back(std::vector<GammaT>{});
	gammav.push_back(std::vector<GammaT>{});
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto t = _Cmplx(0.10L * i, 0.10L * j);
	    zv.back().push_back(t);
	    gammav.back().push_back(lanczos_log_gamma1p(t - GammaT{1}));
	  }
      }

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto z = zv[i - i_min][j - j_min];
	    auto gamma = gammav[i - i_min][j - j_min];
	    data << std::setw(w) << std::real(z)
		 << std::setw(w) << std::imag(z)
		 << std::setw(w) << std::real(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto z = zv[i - i_min][j - j_min];
	    auto gamma = gammav[i - i_min][j - j_min];
	    data << std::setw(w) << std::real(z)
		 << std::setw(w) << std::imag(z)
		 << std::setw(w) << std::imag(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto z = zv[i - i_min][j - j_min];
	    auto gamma = gammav[i - i_min][j - j_min];
	    data << std::setw(w) << std::real(z)
		 << std::setw(w) << std::imag(z)
		 << std::setw(w) << std::abs(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto z = zv[i - i_min][j - j_min];
	    auto gamma = gammav[i - i_min][j - j_min];
	    data << std::setw(w) << std::real(z)
		 << std::setw(w) << std::imag(z)
		 << std::setw(w) << deg * std::arg(gamma)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';
  }

int
main(int n_app_args, char** arg)
{
  std::string plot_data_dir = ".";
  if (n_app_args > 1)
    plot_data_dir = arg[1];

  std::cout << "\n\nLanczos Algorithm\n\n";
  std::cout << "\nlanczos<float>\n";
  plot_lanczos<float>(plot_data_dir + '/' + "log_gamma_lanczos_float.txt");
  std::cout << "\nlanczos<double>\n";
  plot_lanczos<double>(plot_data_dir + '/' + "log_gamma_lanczos_double.txt");
  std::cout << "\nlanczos<long double>\n";
  plot_lanczos<long double>(plot_data_dir + '/' + "log_gamma_lanczos_long_double.txt");
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //std::cout << "\nlanczos<__float128>\n";
  //plot_lanczos<__float128>(plot_data_dir + '/' + "log_gamma_lanczosfloat128.txt");
#endif

  std::cout << "\n\nSpouge Algorithm\n\n";
  std::cout << "\nspouge<float>\n";
  plot_spouge<float>(plot_data_dir + '/' + "log_gamma_spouge_float.txt");
  std::cout << "\nspouge<double>\n";
  plot_spouge<double>(plot_data_dir + '/' + "log_gamma_spouge_double.txt");
  std::cout << "\nspouge<long double>\n";
  plot_spouge<long double>(plot_data_dir + '/' + "log_gamma_spouge_long_double.txt");
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //std::cout << "\nspouge<__float128>\n";
  //plot_spouge<__float128>(plot_data_dir + '/' + "log_gamma_spougefloat128.txt");
#endif
}
