/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <limits>

#include <emsr/float128_io.h>
#include <emsr/float128_math.h>
#include <emsr/special_functions.h>

template<typename Tp>
  void
  test_cyl_bessel()
  {
    using _Val = Tp;
    using _Real = emsr::num_traits_t<_Val>;
    const auto s_pi_2 = emsr::pi_v<_Real> / _Real{2};

    std::cout.precision(emsr::digits10<_Real>());
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    emsr::cyl_neumann(1.0, 0.01);

    std::vector<Tp> nuvec
    {
      Tp{0},
      Tp{1} / Tp{3},
      Tp{1} / Tp{2},
      Tp{2} / Tp{3},
      Tp{1},
      Tp{2},
      Tp{5},
      Tp{10},
      Tp{20},
      Tp{50},
      Tp{100},
    };

    for (auto nu : nuvec)
      {
	std::cout << "\n\n  nu = " << nu << '\n';
	std::cout << ' ' << std::setw(w) << "x";
	std::cout << ' ' << std::setw(w) << "J_nu(x)";
	std::cout << ' ' << std::setw(w) << "J'_nu(x)";
	std::cout << ' ' << std::setw(w) << "N_nu(x)";
	std::cout << ' ' << std::setw(w) << "N'_nu(x)";
	std::cout << ' ' << std::setw(w) << "pi x W{J,N} / 2";
	std::cout << ' ' << std::setw(w) << "I_nu(x)";
	std::cout << ' ' << std::setw(w) << "I'_nu(x)";
	std::cout << ' ' << std::setw(w) << "K_nu(x)";
	std::cout << ' ' << std::setw(w) << "K'_nu(x)";
	std::cout << ' ' << std::setw(w) << "-x W{I,K}";
	std::cout << '\n';
	const auto del = Tp{1} / Tp{100};
	for (unsigned int i = 0; i <= 1000; ++i)
          {
            auto x = i * del;
            std::cout << ' ' << std::setw(w) << x;
            try
              {
        	auto Cyl = emsr::detail::cyl_bessel_jn(nu, x);
        	std::cout << ' ' << std::setw(w) << Cyl.J_value;
        	std::cout << ' ' << std::setw(w) << Cyl.J_deriv;
        	std::cout << ' ' << std::setw(w) << Cyl.N_value;
        	std::cout << ' ' << std::setw(w) << Cyl.N_deriv;
        	std::cout << ' ' << std::setw(w) << s_pi_2 * x * Cyl.Wronskian();
              }
            catch (std::exception& err)
              {
        	std::cerr << err.what() << '\n';
              }
            try
              {
        	auto Cyl = emsr::detail::cyl_bessel_ik(nu, x);
        	std::cout << ' ' << std::setw(w) << Cyl.I_value;
        	std::cout << ' ' << std::setw(w) << Cyl.I_deriv;
        	std::cout << ' ' << std::setw(w) << Cyl.K_value;
        	std::cout << ' ' << std::setw(w) << Cyl.K_deriv;
        	std::cout << ' ' << std::setw(w) << -x * Cyl.Wronskian();
              }
            catch (std::exception& err)
              {
        	std::cerr << err.what() << '\n';
              }
            std::cout << '\n';
          }
      }

    return;
  }

template<typename Tp>
  void
  test_std_bessel()
  {
    using _Val = Tp;
    using _Real = emsr::num_traits_t<_Val>;

    std::cout.precision(emsr::digits10<_Real>());
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();


    std::vector<Tp> nuvec
    {
      Tp{0},
      Tp{1} / Tp{3},
      Tp{1} / Tp{2},
      Tp{2} / Tp{3},
      Tp{1},
      Tp{2},
      Tp{5},
      Tp{10},
      Tp{20},
      Tp{50},
      Tp{100},
    };

    for (auto nu : nuvec)
      {
	std::cout << "\n\n  nu = " << nu << '\n';
	std::cout << ' ' << std::setw(w) << "x";
	std::cout << ' ' << std::setw(w) << "J_nu(x)";
	std::cout << ' ' << std::setw(w) << "N_nu(x)";
	std::cout << ' ' << std::setw(w) << "I_nu(x)";
	std::cout << ' ' << std::setw(w) << "K_nu(x)";
	const auto del = Tp{1} / Tp{10};
	for (int i = 0; i <= 100; ++i)
	  {
	    auto x = i * del;
	    std::cout << ' ' << std::setw(w) << x;
	    try
	      {
		auto jx = emsr::cyl_bessel_j(nu, x);
		auto nx = emsr::cyl_neumann(nu, x);
		std::cout << ' ' << std::setw(w) << jx;
		std::cout << ' ' << std::setw(w) << nx;
	      }
	    catch (std::exception& err)
	      {
		std::cerr << err.what() << '\n';
		auto ix = emsr::cyl_bessel_i(nu, x);
		auto kx = emsr::cyl_bessel_k(nu, x);
		std::cout << ' ' << std::setw(w) << ix;
		std::cout << ' ' << std::setw(w) << kx;
	      }
	    std::cout << '\n';
	  }
      }
  }


///
///
///
int
main()
{
  test_cyl_bessel<float>();
  test_cyl_bessel<double>();
  test_cyl_bessel<long double>();
  //test_cyl_bessel<__float128>();

  test_std_bessel<float>();
  test_std_bessel<double>();
  test_std_bessel<long double>();
  //test_std_bessel<__float128>();

  return 0;
}

