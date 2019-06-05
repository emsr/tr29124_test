/**
 *
 */

#include <cassert>
#include <ext/cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <complex>
#include <string>
#include <bits/float128_io.h>

#include <wrap_burkhardt.h>
#include <sf_struve.h>

/**
 * Take a hard look at the series/asymptotic transition.
 */
template<typename _Tp>
  void
  test_struve_transition(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto del = _Tp{1} / _Tp{100};
    for (int i = 500; i <= +5500; ++i)
      {
	auto t = del * i;
	std::cout << std::setw(width) << t;
	const auto ndel = _Tp{1};
	for (int n = 0; n <= 5; ++n)
	  {
	    auto nu = ndel * n;
	    auto series = std::__detail::__struve_series<std::__detail::_StruveH>(nu, t);
	    auto asymp = std::__detail::__struve_asymp<std::__detail::_StruveK>(nu, t)
		       + std::__detail::__cyl_neumann_n(nu, t);
	    std::cout << '\t'
		      << std::setw(width) << series
		      << std::setw(width) << asymp
		      << std::setw(width) << asymp - series;
	  }
	std::cout << '\n';
      }
    std::cout << "\n\n";
  }


/**
 * Plot the Struve functions.
 */
template<typename _Tp>
  void
  plot_struve(std::string filename, _Tp proto = _Tp{})
  {
    auto data = std::ofstream(filename);

    std::cout.precision(__gnu_cxx::__digits10(proto));
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    data << "\n\n";
    data << "#"
	 << std::setw(width) << "t"
	 << std::setw(width) << "H"
	 << std::setw(width) << "L"
	 << '\n';
    data << "#"
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << '\n';
    const auto del = _Tp{1} / _Tp{100};
    for (int i = 0; i <= +3000; ++i)
      {
	auto t = del * i;
	data << std::setw(width) << t;
	const auto ndel = _Tp{1};
	for (int n = 0; n <= 20; ++n)
	  {
	    auto nu = ndel * n;
	    data << '\t'
		 << std::setw(width) << __gnu_cxx::struve_h(nu, t)
		 << std::setw(width) << __gnu_cxx::struve_l(nu, t);
	  }
	data << '\n';
      }

    data << "\n\n";
    data << "#"
	 << std::setw(width) << "t"
	 << std::setw(width) << "K"
	 << std::setw(width) << "M"
	 << '\n';
    data << "#"
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << std::setw(width) << "========="
	 << '\n';
    for (int i = 0; i <= +3000; ++i)
      {
	auto t = del * i;
	data << std::setw(width) << t;
	const auto ndel = _Tp{1};
	for (int n = 0; n <= 20; ++n)
	  {
	    auto nu = ndel * n;
	    data << '\t'
		 << std::setw(width) << __gnu_cxx::struve_k(nu, t)
		 << std::setw(width) << __gnu_cxx::struve_m(nu, t);
	  }
	data << '\n';
      }
    data << "\n\n";
  }

void
test_struve()
{
  std::cout.precision(__gnu_cxx::__digits10<double>());
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();
  const auto del = double{1} / double{100};
  for (int i = 0; i <= +3000; ++i)
    {
      auto t = del * i;
      std::cout << std::setw(width) << t;
      const auto ndel = double{1};
      for (int n = 0; n <= 20; ++n)
	{
	  auto nu = ndel * n;
	  auto h = __gnu_cxx::struve_h(nu, t);
	  auto l = __gnu_cxx::struve_l(nu, t);
	  auto hb = burkhardt::struve_h(nu, t);
	  auto lb = burkhardt::struve_l(nu, t);
	  std::cout << '\t'
	       << std::setw(width) << h
	       << std::setw(width) << hb
	       << std::setw(width) << (h - hb) / std::abs(hb)
	       << std::setw(width) << l
	       << std::setw(width) << lb
	       << std::setw(width) << (l - lb) / std::abs(lb);
	}
      std::cout << '\n';
    }
  std::cout << "\n\n";
}

int
main(int n_app_args, char** arg)
{
  std::string plot_data_dir = ".";
  if (n_app_args > 1)
    plot_data_dir = arg[1];

  test_struve_transition<float>();
  test_struve_transition<double>();
  test_struve_transition<long double>();
  test_struve_transition<__float128>();

  //using cmplx = std::complex<double>;
  plot_struve<float>(plot_data_dir + '/' + "struve_float.txt");
  plot_struve<double>(plot_data_dir + '/' + "struve_double.txt");
  plot_struve<long double>(plot_data_dir + '/' + "struve_long_double.txt");

  test_struve();

  return 0;
}
