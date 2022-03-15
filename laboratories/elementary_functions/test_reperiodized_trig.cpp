/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <emsr/sf_trig.h>

#include <wrap_boost.h>

template<typename Tp>
  void
  run_sin_cos_pi(Tp proto = Tp{})
  {
    const Tp s_pi = emsr::pi_v<Tp>;

    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sin_pi (GCC)"
	      << std::setw(width) << "sin_pi (Boost)"
	      << std::setw(width) << "delta sin_pi"
	      << std::setw(width) << "delta sin(pi x)"
	      << std::setw(width) << "cos_pi (GCC)"
	      << std::setw(width) << "cos_pi (Boost)"
	      << std::setw(width) << "delta cos_pi"
	      << std::setw(width) << "delta cos(pi x)"
	      << std::setw(width) << "tan_pi (GCC)"
	      << std::setw(width) << "delta tan(pi x)"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    const auto del = Tp{1} / Tp{10};
    for (int i = -1600; i <= +1600; ++i)
      {
	auto x = del * i;
	auto sin_pi_g = emsr::sin_pi(x);
	auto sin_pi_b = beast::sin_pi(x);
	auto cos_pi_g = emsr::cos_pi(x);
	auto cos_pi_b = beast::cos_pi(x);
	auto tan_pi_g = emsr::tan_pi(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << sin_pi_g
		  << std::setw(width) << sin_pi_b
		  << std::setw(width) << sin_pi_g - sin_pi_b
		  << std::setw(width) << sin_pi_g - std::sin(s_pi * x)
		  << std::setw(width) << cos_pi_g
		  << std::setw(width) << cos_pi_b
		  << std::setw(width) << cos_pi_g - cos_pi_b
		  << std::setw(width) << cos_pi_g - std::cos(s_pi * x)
		  << std::setw(width) << tan_pi_g
		  << std::setw(width) << tan_pi_g - std::tan(s_pi * x)
		  << '\n';
      }
    std::cout << '\n';

    std::cout << '\n';
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sin_pi (GCC)"
	      << std::setw(width) << "sin_pi (Boost)"
	      << std::setw(width) << "delta sin_pi"
	      << std::setw(width) << "delta sin(pi x)"
	      << std::setw(width) << "cos_pi (GCC)"
	      << std::setw(width) << "cos_pi (Boost)"
	      << std::setw(width) << "delta cos_pi"
	      << std::setw(width) << "delta cos(pi x)"
	      << std::setw(width) << "tan_pi (GCC)"
	      << std::setw(width) << "delta tan(pi x)"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    const auto bigdel = Tp{1} / Tp{4};
    for (int i = 0; i <= +3200; ++i)
      {
	auto x = bigdel * i;
	auto sin_pi_g = emsr::sin_pi(x);
	auto sin_pi_b = beast::sin_pi(x);
	auto cos_pi_g = emsr::cos_pi(x);
	auto cos_pi_b = beast::cos_pi(x);
	auto tan_pi_g = emsr::tan_pi(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << sin_pi_g
		  << std::setw(width) << sin_pi_b
		  << std::setw(width) << sin_pi_g - sin_pi_b
		  << std::setw(width) << sin_pi_g - std::sin(s_pi * x)
		  << std::setw(width) << cos_pi_g
		  << std::setw(width) << cos_pi_b
		  << std::setw(width) << cos_pi_g - cos_pi_b
		  << std::setw(width) << cos_pi_g - std::cos(s_pi * x)
		  << std::setw(width) << tan_pi_g
		  << std::setw(width) << tan_pi_g - std::tan(s_pi * x)
		  << '\n';
      }
    std::cout << '\n';
  }


int
main()
{
  std::cout << "\ndouble\n=====\n\n";
  run_sin_cos_pi<double>();
}
