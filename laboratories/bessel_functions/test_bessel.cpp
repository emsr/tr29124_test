/*
$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_bessel test_bessel.cpp -lquadmath
./test_bessel > test_bessel.txt

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_bessel test_bessel.cpp -lquadmath
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
    auto w = 8 + std::cout.precision();

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
	const auto del = _Tp{1} / _Tp{100};
	for (unsigned int i = 0; i <= 1000; ++i)
          {
            auto x = i * del;
            std::cout << ' ' << std::setw(w) << x;
            try
              {
        	auto Cyl = std::__detail::__cyl_bessel_jn(nu, x);
        	std::cout << ' ' << std::setw(w) << Cyl.__J_value;
        	std::cout << ' ' << std::setw(w) << Cyl.__J_deriv;
        	std::cout << ' ' << std::setw(w) << Cyl.__N_value;
        	std::cout << ' ' << std::setw(w) << Cyl.__N_deriv;
        	std::cout << ' ' << std::setw(w) << _S_pi_2 * x * Cyl.__Wronskian();
              }
            catch (std::exception& err)
              {
        	std::cerr << err.what() << '\n';
              }
            try
              {
        	auto Cyl = std::__detail::__cyl_bessel_ik(nu, x);
        	std::cout << ' ' << std::setw(w) << Cyl.__I_value;
        	std::cout << ' ' << std::setw(w) << Cyl.__I_deriv;
        	std::cout << ' ' << std::setw(w) << Cyl.__K_value;
        	std::cout << ' ' << std::setw(w) << Cyl.__K_deriv;
        	std::cout << ' ' << std::setw(w) << -x * Cyl.__Wronskian();
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
    auto w = 8 + std::cout.precision();


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
	std::cout << "\n\n  nu = " << nu << '\n';
	std::cout << ' ' << std::setw(w) << "x";
	std::cout << ' ' << std::setw(w) << "J_nu(x)";
	std::cout << ' ' << std::setw(w) << "N_nu(x)";
	std::cout << ' ' << std::setw(w) << "I_nu(x)";
	std::cout << ' ' << std::setw(w) << "K_nu(x)";
	const auto del = _Tp{1} / _Tp{10};
	for (int i = 0; i <= 100; ++i)
	  {
	    auto x = i * del;
	    std::cout << ' ' << std::setw(w) << x;
	    try
	      {
		auto jx = std::cyl_bessel_j(nu, x);
		auto nx = std::cyl_neumann(nu, x);
		std::cout << ' ' << std::setw(w) << jx;
		std::cout << ' ' << std::setw(w) << nx;
	      }
	    catch (std::exception& err)
	      {
		std::cerr << err.what() << '\n';
		auto ix = std::cyl_bessel_i(nu, x);
		auto kx = std::cyl_bessel_k(nu, x);
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
  test_cyl_bessel<__float128>();

  test_std_bessel<float>();
  test_std_bessel<double>();
  test_std_bessel<long double>();
  test_std_bessel<__float128>();

  return 0;
}

