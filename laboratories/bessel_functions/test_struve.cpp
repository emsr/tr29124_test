/**
 *
 */

#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <complex>
#include <string>

#include <emsr/float128_io.h>
#include <emsr/float128_math.h>

#include <wrap_burkhardt.h>
#include <sf_struve.h>

/**
 * Take a hard look at the series/asymptotic transition.
 */
template<typename Tp>
  void
  test_struve_transition(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto del = Tp{1} / Tp{100};
    for (int i = 500; i <= +5500; ++i)
      {
	auto t = del * i;
	std::cout << std::setw(width) << t;
	const auto ndel = Tp{1};
	for (int n = 0; n <= 5; ++n)
	  {
	    auto nu = ndel * n;
	    auto series = emsr::detail::struve_series<emsr::detail::_StruveH>(nu, t);
	    auto asymp = emsr::detail::struve_asymp<emsr::detail::_StruveK>(nu, t)
		       + emsr::detail::cyl_neumann_n(nu, t);
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
template<typename Tp>
  void
  plot_struve(std::string filename, Tp proto = Tp{})
  {
    auto data = std::ofstream(filename);

    std::cout.precision(emsr::digits10(proto));
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
    const auto del = Tp{1} / Tp{100};
    for (int i = 0; i <= +3000; ++i)
      {
	auto t = del * i;
	data << std::setw(width) << t;
	const auto ndel = Tp{1};
	for (int n = 0; n <= 20; ++n)
	  {
	    auto nu = ndel * n;
	    data << '\t'
		 << std::setw(width) << emsr::struve_h(nu, t)
		 << std::setw(width) << emsr::struve_l(nu, t);
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
	const auto ndel = Tp{1};
	for (int n = 0; n <= 20; ++n)
	  {
	    auto nu = ndel * n;
	    data << '\t'
		 << std::setw(width) << emsr::struve_k(nu, t)
		 << std::setw(width) << emsr::struve_m(nu, t);
	  }
	data << '\n';
      }
    data << "\n\n";
  }

void
test_struve()
{
  std::cout.precision(emsr::digits10<double>());
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
	  auto h = emsr::struve_h(nu, t);
	  auto l = emsr::struve_l(nu, t);
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
  //test_struve_transition<__float128>();

  //using cmplx = std::complex<double>;
  plot_struve<float>(plot_data_dir + '/' + "struve_float.txt");
  plot_struve<double>(plot_data_dir + '/' + "struve_double.txt");
  plot_struve<long double>(plot_data_dir + '/' + "struve_long_double.txt");

  test_struve();

  return 0;
}
