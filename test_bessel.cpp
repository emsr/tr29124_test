/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_bessel test_bessel.cpp -lquadmath
./test_bessel > test_bessel.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>
#include <bits/float128_io.h>

template<typename _Tp>
  void
  test_cyl_bessel()
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    const auto _S_pi_2 = __gnu_cxx::__const_pi_half<_Real>();

    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cyl_neumann(1.0, 0.01);

    std::vector<_Tp> nuvec
    {
      _Tp{0},
      _Tp{1} / _Tp{3},
      _Tp{1} / _Tp{2},
      _Tp{2} / _Tp{3},
      _Tp{1},
      _Tp{2},
      _Tp{5},
      _Tp{10},
      _Tp{20},
      _Tp{50},
      _Tp{100},
    };

    for (auto nu : nuvec)
      {
	std::cout << "\n  nu = " << nu << '\n';
	std::cout << ' ' << std::setw(width) << "x";
	std::cout << ' ' << std::setw(width) << "J_nu(x)";
	std::cout << ' ' << std::setw(width) << "J'_nu(x)";
	std::cout << ' ' << std::setw(width) << "N_nu(x)";
	std::cout << ' ' << std::setw(width) << "N'_nu(x)";
	std::cout << ' ' << std::setw(width) << "pi x W{J,N} / 2";
	std::cout << ' ' << std::setw(width) << "I_nu(x)";
	std::cout << ' ' << std::setw(width) << "I'_nu(x)";
	std::cout << ' ' << std::setw(width) << "K_nu(x)";
	std::cout << ' ' << std::setw(width) << "K'_nu(x)";
	std::cout << ' ' << std::setw(width) << "-x W{I,K}";
	std::cout << '\n';
	for (unsigned int i = 0; i <= 1000; ++i)
          {
            auto x = _Tp(i) / _Tp{100};
            std::cout << ' ' << std::setw(width) << x;
            try
              {
        	auto Cyl = std::__detail::__cyl_bessel_jn(nu, x);
        	std::cout << ' ' << std::setw(width) << Cyl.__J_value;
        	std::cout << ' ' << std::setw(width) << Cyl.__J_deriv;
        	std::cout << ' ' << std::setw(width) << Cyl.__N_value;
        	std::cout << ' ' << std::setw(width) << Cyl.__N_deriv;
        	std::cout << ' ' << std::setw(width) << _S_pi_2 * x * Cyl.__Wronskian();
              }
            catch (std::exception& err)
              {
        	std::cerr << err.what() << '\n';
              }
            try
              {
        	auto Cyl = std::__detail::__cyl_bessel_ik(nu, x);
        	std::cout << ' ' << std::setw(width) << Cyl.__I_value;
        	std::cout << ' ' << std::setw(width) << Cyl.__I_deriv;
        	std::cout << ' ' << std::setw(width) << Cyl.__K_value;
        	std::cout << ' ' << std::setw(width) << Cyl.__K_deriv;
        	std::cout << ' ' << std::setw(width) << -x * Cyl.__Wronskian();
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

template<typename _Tp>
  void
  test_std_bessel()
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;

    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n';

    std::vector<_Tp> nuvec
    {
      _Tp{0},
      _Tp{1} / _Tp{3},
      _Tp{1} / _Tp{2},
      _Tp{2} / _Tp{3},
      _Tp{1},
      _Tp{2},
      _Tp{5},
      _Tp{10},
      _Tp{20},
      _Tp{50},
      _Tp{100},
    };

    for (auto nu : nuvec)
      {
	std::cout << ' ' << std::setw(width) << "x";
	std::cout << ' ' << std::setw(width) << "J_nu(x)";
	std::cout << ' ' << std::setw(width) << "N_nu(x)";
	std::cout << ' ' << std::setw(width) << "I_nu(x)";
	std::cout << ' ' << std::setw(width) << "K_nu(x)";
	for (int i = 0; i <= 100; ++i)
	  {
	    auto x = i * _Tp{0.1Q};
	    std::cout << ' ' << std::setw(width) << x;
	    try
	      {
		auto jx = std::cyl_bessel_j(nu, x);
		auto nx = std::cyl_neumann(nu, x);
		std::cout << ' ' << std::setw(width) << jx;
		std::cout << ' ' << std::setw(width) << nx;
	      }
	    catch (std::exception& err)
	      {
		std::cerr << err.what() << '\n';
		auto ix = std::cyl_bessel_i(nu, x);
		auto kx = std::cyl_bessel_k(nu, x);
		std::cout << ' ' << std::setw(width) << ix;
		std::cout << ' ' << std::setw(width) << kx;
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
  test_cyl_bessel<__float128>();

  test_std_bessel<float>();
  test_std_bessel<double>();
  test_std_bessel<long double>();
  test_std_bessel<__float128>();

  return 0;
}

