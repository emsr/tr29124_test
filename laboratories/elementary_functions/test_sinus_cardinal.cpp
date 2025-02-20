/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/float128_io.h>
#include <emsr/sf_cardinal.h>

#include <wrap_boost.h>
#include <wrap_gsl.h>


template<typename Tp>
  void
  test_sinc(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    const auto pi = emsr::pi_v<Tp>;

    std::cout << '\n';
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sinc"
	      << std::setw(width) << "sinc GSL"
	      << std::setw(width) << "sinc Boost"
	      << std::setw(width) << "delta GSL"
	      << std::setw(width) << "delta Boost"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    const auto del = Tp{1} / Tp{10};
    for (int i = -40; i <= +40; ++i)
      {
	auto x = del * i * pi;
	auto sinc = emsr::sinc(x);
	auto sinc_gsl = gsl::sinc(x);
	auto sinc_boost = beast::sinc(x);
	auto delta_gsl = sinc - sinc_gsl;
	auto delta_boost = sinc - sinc_boost;
	std::cout << std::setw(width) << x
		  << std::setw(width) << sinc
		  << std::setw(width) << sinc_gsl
		  << std::setw(width) << sinc_boost
		  << std::setw(width) << delta_gsl
		  << std::setw(width) << delta_boost
		  << '\n';
      }
  }

template<typename Tp>
  void
  test_sinc_pi(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sinc_pi"
	      << std::setw(width) << "sinc_pi GSL"
	      << std::setw(width) << "sinc_pi Boost"
	      << std::setw(width) << "delta GSL"
	      << std::setw(width) << "delta Boost"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    const auto del = Tp{1} / Tp{10};
    for (int i = -40; i <= +40; ++i)
      {
	auto x = del * i;
	auto sinc_pi = emsr::sinc_pi(x);
	auto sinc_pi_gsl = gsl::sinc_pi(x);
	auto sinc_pi_boost = beast::sinc_pi(x);
	auto delta_gsl = sinc_pi - sinc_pi_gsl;
	auto delta_boost = sinc_pi - sinc_pi_boost;
	std::cout << std::setw(width) << x
		  << std::setw(width) << sinc_pi
		  << std::setw(width) << sinc_pi_gsl
		  << std::setw(width) << sinc_pi_boost
		  << std::setw(width) << delta_gsl
		  << std::setw(width) << delta_boost
		  << '\n';
      }
  }

template<typename Tp>
  void
  test_sinhc(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sinhc"
	      << std::setw(width) << "sinhc GSL"
	      << std::setw(width) << "sinhc Boost"
	      << std::setw(width) << "delta GSL"
	      << std::setw(width) << "delta Boost"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';

    const auto del = Tp{1} / Tp{10};
    for (int i = -40; i <= +40; ++i)
      {
	auto x = del * i;
	auto sinhc = emsr::sinhc(x);
	auto sinhc_gsl = gsl::sinhc(x);
	auto sinhc_boost = beast::sinhc(x);
	auto delta_gsl = sinhc - sinhc_gsl;
	auto delta_boost = sinhc - sinhc_boost;
	std::cout << std::setw(width) << x
		  << std::setw(width) << sinhc
		  << std::setw(width) << sinhc_gsl
		  << std::setw(width) << sinhc_boost
		  << std::setw(width) << delta_gsl
		  << std::setw(width) << delta_boost
		  << '\n';
      }
  }

template<typename Tp>
  void
  test_sinhc_pi(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sinhc_pi"
	      << std::setw(width) << "sinhc_pi GSL"
	      << std::setw(width) << "sinhc_pi Boost"
	      << std::setw(width) << "delta GSL"
	      << std::setw(width) << "delta Boost"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    const auto del = Tp{1} / Tp{10};
    for (int i = -40; i <= +40; ++i)
      {
	auto x = del * i;
	auto sinhc_pi = emsr::sinhc_pi(x);
	auto sinhc_pi_gsl = gsl::sinhc_pi(x);
	auto sinhc_pi_boost = beast::sinhc_pi(x);
	auto delta_gsl = sinhc_pi - sinhc_pi_gsl;
	auto delta_boost = sinhc_pi - sinhc_pi_boost;
	std::cout << std::setw(width) << x
		  << std::setw(width) << sinhc_pi
		  << std::setw(width) << sinhc_pi_gsl
		  << std::setw(width) << sinhc_pi_boost
		  << std::setw(width) << delta_gsl
		  << std::setw(width) << delta_boost
		  << '\n';
      }
  }

int
main()
{
  test_sinc<double>();

  test_sinc_pi<double>();

  test_sinhc<double>();

  test_sinhc_pi<double>();
}
