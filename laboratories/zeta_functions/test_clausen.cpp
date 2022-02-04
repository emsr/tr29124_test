/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <emsr/special_functions.h>

#include <wrap_gsl.h>


template<typename Tp>
  void
  test_clausen_cl()
  {
    using cmplx = std::complex<Tp>;

    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << std::setw(width) << "t"
	      << std::setw(width) << "Re(Cl_1)"
	      << std::setw(width) << "Im(Cl_1)"
	      << std::setw(width) << "Re(Cl_2)"
	      << std::setw(width) << "Im(Cl_2)"
	      << std::setw(width) << "Cl_2 GSL"
	      << std::setw(width) << "C_2(x)"
	      << std::setw(width) << "Cl_2"
	      << '\n';
    std::cout << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << '\n';
    const auto del = Tp{1} / Tp{100};
    for (int i = -1000; i <= +1000; ++i)
      {
	auto w = cmplx{del * i};
	auto clausen1 = emsr::clausen(1, w);
	auto clausen2 = emsr::clausen(2, w);
	std::cout << std::setw(width) << std::real(w)
		  << std::setw(width) << std::imag(clausen1)
		  << std::setw(width) << std::real(clausen1)
		  << std::setw(width) << std::real(clausen2)
		  << std::setw(width) << std::imag(clausen2)
		  << std::setw(width) << gsl::clausen_cl(2, std::real(w))
		  << std::setw(width) << emsr::clausen(2, std::real(w))
		  << std::setw(width) << emsr::clausen_cl(2, std::real(w))
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }

template<typename Tp>
  void
  plot_clausen()
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << std::setw(width) << "t"
	      << std::setw(width) << "Cl_1"
	      << std::setw(width) << "Cl_2"
	      << std::setw(width) << "Cl_3"
	      << std::setw(width) << "Cl_4"
	      << std::setw(width) << "Cl_5"
	      << '\n';
    std::cout << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << std::setw(width) << "========="
	      << '\n';
    const auto del = Tp{1} / Tp{100};
    for (int i = -1000; i <= +1000; ++i)
      {
	auto w = del * i;
	std::cout << std::setw(width) << w
		  << std::setw(width) << emsr::clausen_cl(1, w)
		  << std::setw(width) << emsr::clausen_cl(2, w)
		  << std::setw(width) << emsr::clausen_cl(3, w)
		  << std::setw(width) << emsr::clausen_cl(4, w)
		  << std::setw(width) << emsr::clausen_cl(5, w)
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }


int
main()
{
  test_clausen_cl<float>();

  test_clausen_cl<double>();

  test_clausen_cl<long double>();

  //test_clausen_cl<__float128>();

  plot_clausen<double>();
}
