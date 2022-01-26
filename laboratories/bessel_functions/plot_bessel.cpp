/**
 *
 */

#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <cmath>

#include <emsr/special_functions.h>

template<typename _Tp, typename _Bessel>
  void
  plot_bessel(std::string filename, _Bessel bessel)
  {
    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Tp>::digits10);
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    for (int n = 0; n <= 300; ++n)
      {
        auto nu = _Tp(n) / _Tp{3};
	for (int i = 0; i <= 1000; ++i)
	  {
	    auto x = _Tp(0.10L * i);
	    auto j = bessel(nu, x);
	    data << ' ' << std::setw(w) << nu
		 << ' ' << std::setw(w) << x
		 << ' ' << std::setw(w) << j << '\n';
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

  plot_bessel<float>(plot_data_dir + '/' + "cyl_bessel_j_float.txt", std::cyl_bessel_jf);
  plot_bessel<double>(plot_data_dir + '/' + "cyl_bessel_j_double.txt", std::cyl_bessel_j<double, double>);

  plot_bessel<float>(plot_data_dir + '/' + "cyl_neumann_float.txt", std::cyl_neumannf);
  plot_bessel<double>(plot_data_dir + '/' + "cyl_neumann_double.txt", std::cyl_neumann<double, double>);
}
