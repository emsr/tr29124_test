/**
 *
 */

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

#include <emsr/numeric_limits.h>
#include <emsr/sf_bessel.h>
#include <emsr/sf_mod_bessel.h>

template<typename Tp>
  void
  test_cyl_bessel()
  {
    std::cout.precision(emsr::digits10<Tp>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto pi = Tp{3.1415926535897932384626433832795029L};
    auto wronski = [pi](Tp, Tp x) ->Tp { return 2 / (pi * x); };

    std::vector<Tp> nu_vals{120, 128, 136, 160, 200};
    for (auto nu : nu_vals)
      {
	std::cout << '\n';
	for (int i = -10; i <= +10; ++i)
	  {
	    auto x = nu + i * nu / Tp{100};
	    auto bessel = emsr::detail::cyl_bessel_jn(nu, x);
	    auto Wtrue = wronski(nu, x);
	    auto Wcalc = bessel.Wronskian();
	    std::cout << ' ' << std::setw(w) << Wtrue
		      << ' ' << std::setw(w) << Wcalc
		      << ' ' << std::setw(w) << (Wcalc - Wtrue) / Wtrue
		      << ' ' << std::setw(w) << bessel.nu_arg
		      << ' ' << std::setw(w) << bessel.x_arg
		      << ' ' << std::setw(w) << bessel.J_value
		      << ' ' << std::setw(w) << bessel.J_deriv
		      << ' ' << std::setw(w) << bessel.N_value
		      << ' ' << std::setw(w) << bessel.N_deriv
		      << '\n';
	  }
      }
  }

template<typename Tp>
  void
  test_mod_cyl_bessel()
  {
    std::cout.precision(emsr::digits10<Tp>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto wronski = [](Tp, Tp x) ->Tp { return -1 / x; };

    std::vector<Tp> nu_vals{120, 128, 136, 160, 200};
    for (auto nu : nu_vals)
      {
	std::cout << '\n';
	for (int i = -10; i <= +10; ++i)
	  {
	    auto x = nu + i * nu / Tp{100};
	    auto bessel = emsr::detail::cyl_bessel_ik(nu, x);
	    auto Wtrue = wronski(nu, x);
	    auto Wcalc = bessel.Wronskian();
	    std::cout << ' ' << std::setw(w) << Wtrue
		      << ' ' << std::setw(w) << Wcalc
		      << ' ' << std::setw(w) << (Wcalc - Wtrue) / Wtrue
		      << ' ' << std::setw(w) << bessel.nu_arg
		      << ' ' << std::setw(w) << bessel.x_arg
		      << ' ' << std::setw(w) << bessel.I_value
		      << ' ' << std::setw(w) << bessel.I_deriv
		      << ' ' << std::setw(w) << bessel.K_value
		      << ' ' << std::setw(w) << bessel.K_deriv
		      << '\n';
	  }
      }
  }

template<typename Tp>
  void
  test_sph_bessel()
  {
    std::cout.precision(emsr::digits10<Tp>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto wronski = [](int, Tp x) ->Tp { return 1 / (x * x); };

    std::vector<int> n_vals{120, 128, 136, 160, 200};
    for (auto n : n_vals)
      {
	std::cout << '\n';
	for (int i = -10; i <= +10; ++i)
	  {
	    auto x = n + i * n / Tp{100};
	    auto bessel = emsr::detail::sph_bessel_jn(n, x);
	    auto Wtrue = wronski(n, x);
	    auto Wcalc = bessel.Wronskian();
	    std::cout << ' ' << std::setw(w) << Wtrue
		      << ' ' << std::setw(w) << Wcalc
		      << ' ' << std::setw(w) << (Wcalc - Wtrue) / Wtrue
		      << ' ' << std::setw(w) << bessel.n_arg
		      << ' ' << std::setw(w) << bessel.x_arg
		      << ' ' << std::setw(w) << bessel.j_value
		      << ' ' << std::setw(w) << bessel.j_deriv
		      << ' ' << std::setw(w) << bessel.n_value
		      << ' ' << std::setw(w) << bessel.n_deriv
		      << '\n';
	  }
      }
  }

template<typename Tp>
  void
  test_mod_sph_bessel()
  {
    std::cout.precision(emsr::digits10<Tp>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto pi = Tp{3.1415926535897932384626433832795029L};
    auto wronski = [pi](Tp, Tp x) ->Tp { return -pi / (2 * x * x); };

    std::vector<Tp> n_vals{120, 128, 136, 160, 200};
    for (auto n : n_vals)
      {
	std::cout << '\n';
	for (int i = -10; i <= +10; ++i)
	  {
	    auto x = n + i * n / Tp{100};
	    auto bessel = emsr::detail::sph_bessel_ik(n, x);
	    auto Wtrue = wronski(n, x);
	    auto Wcalc = bessel.Wronskian();
	    std::cout << ' ' << std::setw(w) << Wtrue
		      << ' ' << std::setw(w) << Wcalc
		      << ' ' << std::setw(w) << (Wcalc - Wtrue) / Wtrue
		      << ' ' << std::setw(w) << bessel.n_arg
		      << ' ' << std::setw(w) << bessel.x_arg
		      << ' ' << std::setw(w) << bessel.i_value
		      << ' ' << std::setw(w) << bessel.i_deriv
		      << ' ' << std::setw(w) << bessel.k_value
		      << ' ' << std::setw(w) << bessel.k_deriv
		      << '\n';
	  }
      }
  }

int
main()
{
  std::cout << "\n\nCylindrical Bessel functions\n";
  test_cyl_bessel<double>();

  std::cout << "\n\nModified cylindrical Bessel functions\n";
  test_mod_cyl_bessel<double>();

  std::cout << "\n\nSpherical Bessel functions\n";
  test_sph_bessel<double>();

  std::cout << "\n\nModified spherical Bessel functions\n";
  test_mod_sph_bessel<double>();
}
