/*
$HOME/bin/bin/g++ -std=c++17 -I. -o test_large_order_bessel test_large_order_bessel.cpp 
./test_large_order_bessel > test_large_order_bessel.txt
*/

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

template<typename _Tp>
  void
  test_cyl_bessel()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto pi = _Tp{3.1415926535897932384626433832795029L};
    auto wronski = [pi](_Tp nu, _Tp x) ->_Tp { return 2 / (pi * x); };

    std::vector<_Tp> nu_vals{120, 128, 136, 160, 200};
    for (auto nu : nu_vals)
      {
	std::cout << '\n';
	for (int i = -10; i <= +10; ++i)
	  {
	    auto x = nu + i * nu / _Tp{100};
	    auto bessel = std::__detail::__cyl_bessel_jn(nu, x);
	    auto Wtrue = wronski(nu, x);
	    auto Wcalc = bessel.__Wronskian();
	    std::cout << ' ' << std::setw(w) << Wtrue
		      << ' ' << std::setw(w) << Wcalc
		      << ' ' << std::setw(w) << (Wcalc - Wtrue) / Wtrue
		      << ' ' << std::setw(w) << bessel.__nu_arg
		      << ' ' << std::setw(w) << bessel.__x_arg
		      << ' ' << std::setw(w) << bessel.__J_value
		      << ' ' << std::setw(w) << bessel.__J_deriv
		      << ' ' << std::setw(w) << bessel.__N_value
		      << ' ' << std::setw(w) << bessel.__N_deriv
		      << '\n';
	  }
      }
  }

template<typename _Tp>
  void
  test_mod_cyl_bessel()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto wronski = [](_Tp, _Tp x) ->_Tp { return -1 / x; };

    std::vector<_Tp> nu_vals{120, 128, 136, 160, 200};
    for (auto nu : nu_vals)
      {
	std::cout << '\n';
	for (int i = -10; i <= +10; ++i)
	  {
	    auto x = nu + i * nu / _Tp{100};
	    auto bessel = std::__detail::__cyl_bessel_ik(nu, x);
	    auto Wtrue = wronski(nu, x);
	    auto Wcalc = bessel.__Wronskian();
	    std::cout << ' ' << std::setw(w) << Wtrue
		      << ' ' << std::setw(w) << Wcalc
		      << ' ' << std::setw(w) << (Wcalc - Wtrue) / Wtrue
		      << ' ' << std::setw(w) << bessel.__nu_arg
		      << ' ' << std::setw(w) << bessel.__x_arg
		      << ' ' << std::setw(w) << bessel.__I_value
		      << ' ' << std::setw(w) << bessel.__I_deriv
		      << ' ' << std::setw(w) << bessel.__K_value
		      << ' ' << std::setw(w) << bessel.__K_deriv
		      << '\n';
	  }
      }
  }

template<typename _Tp>
  void
  test_sph_bessel()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto wronski = [](int, _Tp x) ->_Tp { return 1 / (x * x); };

    std::vector<int> n_vals{120, 128, 136, 160, 200};
    for (auto n : n_vals)
      {
	std::cout << '\n';
	for (int i = -10; i <= +10; ++i)
	  {
	    auto x = n + i * n / _Tp{100};
	    auto bessel = std::__detail::__sph_bessel_jn(n, x);
	    auto Wtrue = wronski(n, x);
	    auto Wcalc = bessel.__Wronskian();
	    std::cout << ' ' << std::setw(w) << Wtrue
		      << ' ' << std::setw(w) << Wcalc
		      << ' ' << std::setw(w) << (Wcalc - Wtrue) / Wtrue
		      << ' ' << std::setw(w) << bessel.__n_arg
		      << ' ' << std::setw(w) << bessel.__x_arg
		      << ' ' << std::setw(w) << bessel.__j_value
		      << ' ' << std::setw(w) << bessel.__j_deriv
		      << ' ' << std::setw(w) << bessel.__n_value
		      << ' ' << std::setw(w) << bessel.__n_deriv
		      << '\n';
	  }
      }
  }

template<typename _Tp>
  void
  test_mod_sph_bessel()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto pi = _Tp{3.1415926535897932384626433832795029L};
    auto wronski = [pi](_Tp, _Tp x) ->_Tp { return -pi / (2 * x * x); };

    std::vector<_Tp> n_vals{120, 128, 136, 160, 200};
    for (auto n : n_vals)
      {
	std::cout << '\n';
	for (int i = -10; i <= +10; ++i)
	  {
	    auto x = n + i * n / _Tp{100};
	    auto bessel = std::__detail::__sph_bessel_ik(n, x);
	    auto Wtrue = wronski(n, x);
	    auto Wcalc = bessel.__Wronskian();
	    std::cout << ' ' << std::setw(w) << Wtrue
		      << ' ' << std::setw(w) << Wcalc
		      << ' ' << std::setw(w) << (Wcalc - Wtrue) / Wtrue
		      << ' ' << std::setw(w) << bessel.__n_arg
		      << ' ' << std::setw(w) << bessel.__x_arg
		      << ' ' << std::setw(w) << bessel.__i_value
		      << ' ' << std::setw(w) << bessel.__i_deriv
		      << ' ' << std::setw(w) << bessel.__k_value
		      << ' ' << std::setw(w) << bessel.__k_deriv
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
