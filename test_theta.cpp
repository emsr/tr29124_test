/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_theta test_theta.cpp -lquadmath
./test_theta > test_theta.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <bits/float128_io.h>

template<typename _Tp>
  void
  test_theta(_Tp proto = _Tp{})
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto width = 8 + std::cout.precision();
    std::cout << std::showpoint << std::scientific;

    std::cout << "\n\n Theta function values\n";
    std::cout << " =====================\n";
    for (int i = 0; i <= 20; ++i)
      {
	auto nu = _Tp(i * 0.1Q);
	std::cout << '\n' << " nu   = " << std::setw(width) << nu << '\n';
	std::cout << ' ' << std::setw(width) << "x"
		  << ' ' << std::setw(width) << "theta_1"
		  << ' ' << std::setw(width) << "theta_2"
		  << ' ' << std::setw(width) << "theta_3"
		  << ' ' << std::setw(width) << "theta_4"
		  << '\n';
	std::cout << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << '\n';
	for (int j = 0; j <= 100; ++j)
	  {
	    auto x = _Tp(j * 0.01Q);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << __gnu_cxx::theta_1(nu, x)
		      << ' ' << std::setw(width) << __gnu_cxx::theta_2(nu, x)
		      << ' ' << std::setw(width) << __gnu_cxx::theta_3(nu, x)
		      << ' ' << std::setw(width) << __gnu_cxx::theta_4(nu, x)
		      << '\n';
	  }
      }

    std::cout << "\n\n Theta function values and compares with Jacobi elliptic functions\n";
    std::cout << " =================================================================\n";
    auto k = _Tp{1} / _Tp{3};
    std::cout << '\n' << " k    = " << std::setw(width) << k;
    std::cout << '\n' << " q(k) = " << std::setw(width) << __gnu_cxx::ellnome(k) << '\n';
    std::cout << ' ' << std::setw(width) << "x"
	      << ' ' << std::setw(width) << "theta_s"
	      << ' ' << std::setw(width) << "theta_c"
	      << ' ' << std::setw(width) << "theta_n"
	      << ' ' << std::setw(width) << "theta_d"
	      << ' ' << std::setw(width) << "th_s/th_n - sn"
	      << ' ' << std::setw(width) << "th_c/th_n - cn"
	      << ' ' << std::setw(width) << "th_d/th_n - dn"
	      << '\n';
    std::cout << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << ' ' << std::setw(width) << "------------------------"
	      << '\n';
    for (int j = -1000; j <= 1000; ++j)
      {
	auto x = _Tp(j * 0.01Q);
	auto s = __gnu_cxx::theta_s(k, x);
	auto c = __gnu_cxx::theta_c(k, x);
	auto d = __gnu_cxx::theta_d(k, x);
	auto n = __gnu_cxx::theta_n(k, x);
	auto sncndn = std::__detail::__jacobi_sncndn(k, x);
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << s
		  << ' ' << std::setw(width) << c
		  << ' ' << std::setw(width) << d
		  << ' ' << std::setw(width) << n
		  << ' ' << std::setw(width) << s / n - std::get<0>(sncndn)
		  << ' ' << std::setw(width) << c / n - std::get<1>(sncndn)
		  << ' ' << std::setw(width) << d / n - std::get<2>(sncndn)
		  << '\n';
      }

    std::cout << "\n\n Theta function values and compares with Jacobi elliptic functions\n";
    std::cout << " =================================================================\n";
    for (int i = -10; i <= 10; ++i)
      {
	auto k = _Tp(i * 0.1Q);
	std::cout << '\n' << " k    = " << std::setw(width) << k;
	std::cout << '\n' << " q(k) = " << std::setw(width) << __gnu_cxx::ellnome(k) << '\n';
	std::cout << ' ' << std::setw(width) << "x"
		  << ' ' << std::setw(width) << "theta_s"
		  << ' ' << std::setw(width) << "theta_c"
		  << ' ' << std::setw(width) << "theta_n"
		  << ' ' << std::setw(width) << "theta_d"
		  << ' ' << std::setw(width) << "th_s/th_n - sn"
		  << ' ' << std::setw(width) << "th_c/th_n - cn"
		  << ' ' << std::setw(width) << "th_d/th_n - dn"
		  << '\n';
	std::cout << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << '\n';
	for (int j = 0; j <= 100; ++j)
	  {
	    auto x = _Tp(j * 0.01Q);
	    auto s = __gnu_cxx::theta_s(k, x);
	    auto c = __gnu_cxx::theta_c(k, x);
	    auto d = __gnu_cxx::theta_d(k, x);
	    auto n = __gnu_cxx::theta_n(k, x);
	    auto sncndn = std::__detail::__jacobi_sncndn(k, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << s
		      << ' ' << std::setw(width) << c
		      << ' ' << std::setw(width) << d
		      << ' ' << std::setw(width) << n
		      << ' ' << std::setw(width) << s / n - std::get<0>(sncndn)
		      << ' ' << std::setw(width) << c / n - std::get<1>(sncndn)
		      << ' ' << std::setw(width) << d / n - std::get<2>(sncndn)
		      << '\n';
	  }
      }
  }

int
main()
{
  std::cout << "\nfloat\n===========\n";
  test_theta<float>();
  std::cout << "\ndouble\n===========\n";
  test_theta<double>();
  std::cout << "\nlong double\n===========\n";
  test_theta<long double>();
  std::cout << "\n__float128\n===========\n";
  test_theta<__float128>();
}
