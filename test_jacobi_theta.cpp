/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_theta test_jacobi_theta.cpp -lquadmath
./test_jacobi_theta > test_jacobi_theta.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_theta test_jacobi_theta.cpp -lquadmath
./test_jacobi_theta > test_jacobi_theta.txt
*/

#include <iostream>
#include <iomanip>
#include <ext/cmath>
#include <bits/float128_io.h>

template<typename _Tp>
  void
  test_jacobi_theta(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    std::cout << "\n\n Theta function values\n";
    std::cout << " =====================\n";
    const auto del1 = _Tp{1} / _Tp{10};
    const auto del01 = _Tp{1} / _Tp{100};
    for (int i = 0; i <= 20; ++i)
      {
	auto nu = i * del1;
	std::cout << '\n' << " nu   = " << std::setw(width) << nu << '\n';
	std::cout << ' ' << std::setw(width) << "x"
		  << ' ' << std::setw(width) << "jacobi_theta_1"
		  << ' ' << std::setw(width) << "jacobi_theta_2"
		  << ' ' << std::setw(width) << "jacobi_theta_3"
		  << ' ' << std::setw(width) << "jacobi_theta_4"
		  << '\n';
	std::cout << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << ' ' << std::setw(width) << "------------------------"
		  << '\n';
	for (int j = 0; j <= 100; ++j)
	  {
	    auto x = j * del01;
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << __gnu_cxx::jacobi_theta_1(nu, x)
		      << ' ' << std::setw(width) << __gnu_cxx::jacobi_theta_2(nu, x)
		      << ' ' << std::setw(width) << __gnu_cxx::jacobi_theta_3(nu, x)
		      << ' ' << std::setw(width) << __gnu_cxx::jacobi_theta_4(nu, x)
		      << '\n';
	  }
      }
  }

int
main()
{
  auto jt1 [[maybe_unused]] = std::__detail::__jacobi_theta_1_sum(2.0, -8.0);
  auto jt2 [[maybe_unused]] = std::__detail::__jacobi_theta_2_sum(1.999, -8.0);
  auto jt3 [[maybe_unused]] = std::__detail::__jacobi_theta_3_sum(3.0, -8.0);
  auto jt4 [[maybe_unused]] = std::__detail::__jacobi_theta_4_sum(2.999, -8.0);
  test_jacobi_theta(1.0);
}
