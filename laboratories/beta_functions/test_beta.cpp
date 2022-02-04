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

#include <emsr/float128_io.h>
#include <emsr/numeric_limits.h>
#include <emsr/sf_beta.h>

#include <wrap_boost.h>

template<typename Tp>
  void
  test_beta(Tp proto = Tp{})
  {
    //using _Val = Tp;
    //using _Real = emsr::num_traits_t<_Val>;

    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "b"
	      << ' ' << std::setw(width) << "beta"
	      << ' ' << std::setw(width) << "boost::beta"
	      << ' ' << std::setw(width) << "delta_boost"
	      << '\n';
    int i_min = 1;
    for (int i = i_min; i <= +500; ++i)
      {
	auto a = Tp{0.1Q} * i;
	int j_min = 1;
	for (int j = j_min; j <= +500; ++j)
	  {
	    auto b = Tp{0.1Q} * j;
	    auto gbet = emsr::detail::beta(a, b);
	    auto bbet = beast::beta(a, b);
	    std::cout << ' ' << std::setw(width) << a
		      << ' ' << std::setw(width) << b
		      << ' ' << std::setw(width) << gbet
		      << ' ' << std::setw(width) << bbet
		      << ' ' << std::setw(width) << (gbet - bbet) / std::abs(bbet)
		      << '\n';
	  }
      }
  }

template<typename Tp>
  void
  plot_beta(std::string filename)
  {
    using _Val = Tp;
    using _Real = emsr::num_traits_t<_Val>;

    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Real>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    int i_min = -150;
    for (int i = i_min; i <= +150; ++i)
      {
	auto a = Tp{0.02L} * i;
	int j_min = -150;
	data << '\n';
	for (int j = j_min; j <= +150; ++j)
	  {
	    auto b = Tp{0.02L} * j;
	    auto gbet = emsr::detail::beta(a, b);
	    data << ' ' << std::setw(width) << a
		 << ' ' << std::setw(width) << b
		 << ' ' << std::setw(width) << gbet
		 << '\n';
	  }
      }
  }

int
main(int n_app_args, char** arg)
{
  std::string plot_data_dir = ".";
  if (n_app_args > 1)
    plot_data_dir = arg[1];

  emsr::detail::beta(0.1F, 1.9F);
  emsr::detail::beta(0.1F, 35.1F);

  test_beta<float>();

  test_beta<double>();

  test_beta<long double>();

  //test_beta<__float128>();

  // Beta seems to be either really tiny or really huge.
  // Maybe graph log_beta.
  plot_beta<double>(plot_data_dir + '/' + "beta_double.txt");
}
