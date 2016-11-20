/*
$HOME/bin_tr29124/bin/g++ -std=c++17 -g -Wall -Wextra -Wno-psabi -I. -o plot_bessel plot_bessel.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./plot_bessel

$HOME/bin/bin/g++ -std=gnu++17 -DNO_LOGBQ -g -Wall -Wextra -I. -o plot_bessel plot_bessel.cpp -lquadmath
./plot_bessel > plot_bessel.txt
*/

#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <ext/cmath>
#include <bits/specfun.h>

template<typename _Tp, typename _Bessel>
  void
  plot_bessel(std::string filename, _Bessel __bessel)
  {
    auto data = std::ofstream(filename);

    data.precision(std::numeric_limits<_Tp>::digits10);
    data << std::showpoint << std::scientific;
    auto width = 8 + data.precision();

    for (int n = 0; n <= 300; ++n)
      {
        auto nu = _Tp(n) / _Tp{3};
	for (int i = 0; i <= 1000; ++i)
	  {
	    auto x = _Tp(0.10L * i);
	    auto j = __bessel(nu, x);
	    data << ' ' << std::setw(width) << nu
		 << ' ' << std::setw(width) << x
		 << ' ' << std::setw(width) << j << '\n';
	  }
	data << '\n';
      }
    data << '\n';
  }

int
main()
{
  plot_bessel<float>("plot/cyl_bessel_j_float.txt", std::cyl_bessel_jf);
  plot_bessel<double>("plot/cyl_bessel_j_double.txt", std::cyl_bessel_j<double, double>);

  plot_bessel<float>("plot/cyl_neumann_float.txt", std::cyl_neumannf);
  plot_bessel<double>("plot/cyl_neumann_double.txt", std::cyl_neumann<double, double>);
}
