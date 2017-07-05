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
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    std::cout << "\n\n Theta function values\n";
    std::cout << " =====================\n";
    const auto del1 = _Tp{1} / _Tp{10};
    const auto del01 = _Tp{1} / _Tp{100};
    for (int i = -9; i <= 9; ++i)
      {
	auto q = i * del1;
	std::cout << '\n' << " q   = " << std::setw(w) << q << '\n';
	std::cout << ' ' << std::setw(w) << "x"
		  << ' ' << std::setw(w) << "jacobi_theta_1"
		  << ' ' << std::setw(w) << "jacobi_theta_2"
		  << ' ' << std::setw(w) << "jacobi_theta_3"
		  << ' ' << std::setw(w) << "jacobi_theta_4"
		  << '\n';
	std::cout << ' ' << std::setw(w) << "------------------------"
		  << ' ' << std::setw(w) << "------------------------"
		  << ' ' << std::setw(w) << "------------------------"
		  << ' ' << std::setw(w) << "------------------------"
		  << ' ' << std::setw(w) << "------------------------"
		  << '\n';
	for (int j = 0; j <= 100; ++j)
	  {
	    auto x = j * del01;
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << std::__detail::__jacobi_theta_1(q, x)
		      << ' ' << std::setw(w) << std::__detail::__jacobi_theta_2(q, x)
		      << ' ' << std::setw(w) << std::__detail::__jacobi_theta_3(q, x)
		      << ' ' << std::setw(w) << std::__detail::__jacobi_theta_4(q, x)
		      << '\n';
	  }
      }
  }

template<typename _Tp>
  void
  plot_jacobi_theta(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto _S_pi = __gnu_cxx::__const_pi(proto);

    std::cout << "\n\n\n";
    _Tp q0 = _Tp{0.15L};
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _S_pi * i / 100;
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << std::__detail::__jacobi_theta_1(q0, x)
		  << ' ' << std::setw(w) << std::__detail::__jacobi_theta_2(q0, x)
		  << ' ' << std::setw(w) << std::__detail::__jacobi_theta_3(q0, x)
		  << ' ' << std::setw(w) << std::__detail::__jacobi_theta_4(q0, x)
		  << '\n';
      }

    std::cout << "\n\n\n";
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _S_pi * i / 100;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {_Tp{0.05L}, _Tp{0.5L}, _Tp{0.7L}, _Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_1(q, x);
	std::cout << '\n';
      }

    std::cout << "\n\n\n";
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _S_pi * i / 100;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {_Tp{0.05L}, _Tp{0.5L}, _Tp{0.7L}, _Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_2(q, x);
	std::cout << '\n';
      }

    std::cout << "\n\n\n";
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _S_pi * i / 100;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {_Tp{0.05L}, _Tp{0.5L}, _Tp{0.7L}, _Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_3(q, x);
	std::cout << '\n';
      }

    std::cout << "\n\n\n";
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _S_pi * i / 100;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {_Tp{0.05L}, _Tp{0.5L}, _Tp{0.7L}, _Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_4(q, x);
	std::cout << '\n';
      }

    std::cout << "\n\n\n";
    for (int i = 0; i < 100; ++i)
      {
	auto q = _Tp(i) / _Tp{100.0L};
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {_Tp{0.0L}, _Tp{0.4L}, _Tp{5.0L}, _Tp{10.0L}, _Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_1(q, x);
	std::cout << '\n';
      }

    std::cout << "\n\n\n";
    for (int i = 0; i < 100; ++i)
      {
	auto q = _Tp(i) / _Tp{100.0L};
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {_Tp{0.0L}, _Tp{0.4L}, _Tp{5.0L}, _Tp{10.0L}, _Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_2(q, x);
	std::cout << '\n';
      }

    std::cout << "\n\n\n";
    for (int i = 0; i < 100; ++i)
      {
	auto q = _Tp(i) / _Tp{100.0L};
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {_Tp{0.0L}, _Tp{0.4L}, _Tp{5.0L}, _Tp{10.0L}, _Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_3(q, x);
	std::cout << '\n';
      }

    std::cout << "\n\n\n";
    for (int i = 0; i < 100; ++i)
      {
	auto q = _Tp(i) / _Tp{100.0L};
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {_Tp{0.0L}, _Tp{0.4L}, _Tp{5.0L}, _Tp{10.0L}, _Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_4(q, x);
	std::cout << '\n';
      }
  }


int
main()
{
  plot_jacobi_theta(1.0);

  //auto jt1 [[maybe_unused]] = std::__detail::__jacobi_theta_1_sum(2.0, -8.0);
  //auto jt2 [[maybe_unused]] = std::__detail::__jacobi_theta_2_sum(1.999, -8.0);
  //auto jt3 [[maybe_unused]] = std::__detail::__jacobi_theta_3_sum(3.0, -8.0);
  //auto jt4 [[maybe_unused]] = std::__detail::__jacobi_theta_4_sum(2.999, -8.0);

  test_jacobi_theta(1.0);
}
