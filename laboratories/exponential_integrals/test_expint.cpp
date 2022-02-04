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
#include <emsr/fp_type_util.h>
#include <emsr/numeric_limits.h>
#include <emsr/sf_expint.h>

/*
                 series     asymp     large-n
float                                     200
double                       2          20000
long double                             25000
__float128                                   

large-n works brilliantly for n=0?!?
*/

#include <wrap_boost.h>

template<typename Tp>
  void
  test_expint()
  {
    using _Val = Tp;
    using _Real = emsr::num_traits_t<_Val>;
    auto _S_NaN = emsr::quiet_NaN<_Real>();
    std::cout.precision(emsr::digits10<_Real>());
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
	const auto del = Tp{1} / Tp{10};
	for (int i = i_min; i <= +500; ++i)
	  {
	    auto x = del * i;

	    Tp ens = _S_NaN;
	    try
	    {
	      ens = emsr::detail::expint_En_series(n, x);
	    }
	    catch (...)
	    {
	    }

	    Tp enc = _S_NaN;
	    try
	    {
	      enc = emsr::detail::expint_En_cont_frac(n, x);
	    }
	    catch (...)
	    {
	    }

	    Tp enn = _S_NaN;
	    try
	    {
	      enn = emsr::detail::expint_En_large_n(n, x);
	    }
	    catch (...)
	    {
	    }

	    Tp ena = _S_NaN;
	    try
	    {
	      ena = emsr::detail::expint_En_asymp(n, x);
	    }
	    catch (...)
	    {
	    }

	    Tp enb = _S_NaN;
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

  //test_expint<__float128>();
}
