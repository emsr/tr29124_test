/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_expint test_expint.cpp -lquadmath -Lwrappers/debug -lwrap_boost
./test_expint > test_expint.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_expint test_expint.cpp -lquadmath -Lwrappers/debug -lwrap_boost
./test_expint > test_expint.txt

g++ -std=gnu++17 -DNO_LOGBQ -g -Wall -Wextra -I. -o test_expint test_expint.cpp -lquadmath -Lwrappers -lwrap_boost
./test_expint > test_expint.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>
#include <bits/float128_io.h>

/*
                 series     asymp     large-n
float                                     200
double                       2          20000
long double                             25000
__float128                                   

large-n works brilliantly for n=0?!?
*/

#include "wrap_boost.h"

template<typename _Tp>
  void
  test_expint()
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    auto _S_NaN = __gnu_cxx::__quiet_NaN<_Real>();
    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();
    std::vector<unsigned int>
      order({0, 1, 2, 5, 10, 20, 50, 100, 200, 500,
	  1000, 2000, 5000, 10000, 20000, 50000, 100000});

    for (auto n : order)
      {
	std::cout << '\n'
		  << ' ' << std::setw(4) << "n"
		  << ' ' << std::setw(width) << "x"
		  << ' ' << std::setw(width) << "E_n boost"
		  << ' ' << std::setw(width) << "E_n series"
		  << ' ' << std::setw(width) << "delta series"
		  << ' ' << std::setw(width) << "E_n cfrac"
		  << ' ' << std::setw(width) << "delta cfrac"
		  << ' ' << std::setw(width) << "E_n large n"
		  << ' ' << std::setw(width) << "delta large n"
		  << ' ' << std::setw(width) << "E_n asymp"
		  << ' ' << std::setw(width) << "delta asymp"
		  << '\n';
	int i_min = -500;
	const auto del = _Tp{1} / _Tp{10};
	for (int i = i_min; i <= +500; ++i)
	  {
	    auto x = del * i;

	    _Tp ens = _S_NaN;
	    try
	    {
	      ens = std::__detail::__expint_En_series(n, x);
	    }
	    catch (...)
	    {
	    }

	    _Tp enc = _S_NaN;
	    try
	    {
	      enc = std::__detail::__expint_En_cont_frac(n, x);
	    }
	    catch (...)
	    {
	    }

	    _Tp enn = _S_NaN;
	    try
	    {
	      enn = std::__detail::__expint_En_large_n(n, x);
	    }
	    catch (...)
	    {
	    }

	    _Tp ena = _S_NaN;
	    try
	    {
	      ena = std::__detail::__expint_En_asymp(n, x);
	    }
	    catch (...)
	    {
	    }

	    _Tp enb = _S_NaN;
	    try
	    {
	      enb = beast::expint(n, x);
	    }
	    catch (...)
	    {
	    }

	    std::cout << ' ' << std::setw(4) << n
		      << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << enb
		      << ' ' << std::setw(width) << ens
		      << ' ' << std::setw(width) << (ens - enb) / std::abs(enb)
		      << ' ' << std::setw(width) << enc
		      << ' ' << std::setw(width) << (enc - enb) / std::abs(enb)
		      << ' ' << std::setw(width) << enn
		      << ' ' << std::setw(width) << (enn - enb) / std::abs(enb)
		      << ' ' << std::setw(width) << ena
		      << ' ' << std::setw(width) << (ena - enb) / std::abs(enb)
		      << '\n';
	  }
      }
  }

int
main()
{
  test_expint<float>();

  test_expint<double>();

  test_expint<long double>();

  test_expint<__float128>();
}
