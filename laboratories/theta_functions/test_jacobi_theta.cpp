/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_theta test_jacobi_theta.cpp -lquadmath
./test_jacobi_theta > test_jacobi_theta.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_theta test_jacobi_theta.cpp -lquadmath
./test_jacobi_theta > test_jacobi_theta.txt
*/

#include <iostream>
#include <iomanip>
#include <ext/cmath>
#include <bits/float128_io.h>

/**
 * Plot the thetas over x= [0, 2pi] for q = 0, 0.1, ..., 0.9
 */
template<typename _Tp>
  void
  test_jacobi_theta(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto _S_pi = __gnu_cxx::__const_pi(proto);

    const auto del1 = _Tp{1} / _Tp{10};
    const auto del01 = _Tp{1} / _Tp{100};
    for (int i = 0; i <= 9; ++i)
      {
	auto q = i * del1;
	std::cout << '\n' << '\n' << " q   = " << std::setw(w) << q << '\n';
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
	for (int j = 0; j <= 200; ++j)
	  {
	    auto x = j * del01;
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << std::__detail::__jacobi_theta_1(q, _S_pi * x)
		      << ' ' << std::setw(w) << std::__detail::__jacobi_theta_2(q, _S_pi * x)
		      << ' ' << std::setw(w) << std::__detail::__jacobi_theta_3(q, _S_pi * x)
		      << ' ' << std::setw(w) << std::__detail::__jacobi_theta_4(q, _S_pi * x)
		      << '\n';
	  }
      }
    std::cout.flush();
  }

/**
 * Reproduce DLMF plots.
 */
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
    std::cout << "q = " << q0 << '\n';
    for (int i = 0; i <= 200; ++i)
      {
	auto x = _Tp(i) / _Tp{100};
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << std::__detail::__jacobi_theta_1(q0, _S_pi * x)
		  << ' ' << std::setw(w) << std::__detail::__jacobi_theta_2(q0, _S_pi * x)
		  << ' ' << std::setw(w) << std::__detail::__jacobi_theta_3(q0, _S_pi * x)
		  << ' ' << std::setw(w) << std::__detail::__jacobi_theta_4(q0, _S_pi * x)
		  << '\n';
      }
    std::cout.flush();

    // Vary x for fixed q.
    int nx = 400;
    auto delx = 2 / _Tp(nx);

    std::cout << "\n\n\n";
    std::cout << "theta_1; q = 0.05, 0.5, 0.7, 0.9\n";
    for (int i = 0; i <= nx; ++i)
      {
	auto x = _Tp(i) * delx;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {_Tp{0.05L}, _Tp{0.5L}, _Tp{0.7L}, _Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_1(q, _S_pi * x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_2; q = 0.05, 0.5, 0.7, 0.9\n";
    for (int i = 0; i <= nx; ++i)
      {
	auto x = _Tp(i) * delx;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {_Tp{0.05L}, _Tp{0.5L}, _Tp{0.7L}, _Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_2(q, _S_pi * x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_3; q = 0.05, 0.5, 0.7, 0.9\n";
    for (int i = 0; i <= nx; ++i)
      {
	auto x = _Tp(i) * delx;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {_Tp{0.05L}, _Tp{0.5L}, _Tp{0.7L}, _Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_3(q, _S_pi * x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_4; q = 0.05, 0.5, 0.7, 0.9\n";
    for (int i = 0; i <= nx; ++i)
      {
	auto x = _Tp(i) * delx;
	std::cout << ' ' << std::setw(w) << x;
	for (auto q : {_Tp{0.05L}, _Tp{0.5L}, _Tp{0.7L}, _Tp{0.9L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_4(q, _S_pi * x);
	std::cout << '\n';
      }
    std::cout.flush();

    // Vary q for fixed x.
    int nq = 200;
    auto delq = 1 / _Tp(nq);

    std::cout << "\n\n\n";
    std::cout << "theta_1; x = 0, 0.4, 5, 10, 40\n";
    for (int i = 0; i < nq; ++i)
      {
	auto q = i * delq;
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {_Tp{0.0L}, _Tp{0.4L}, _Tp{5.0L}, _Tp{10.0L}, _Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_1(q, x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_2; x = 0, 0.4, 5, 10, 40\n";
    for (int i = 0; i < nq; ++i)
      {
	auto q = i * delq;
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {_Tp{0.0L}, _Tp{0.4L}, _Tp{5.0L}, _Tp{10.0L}, _Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_2(q, x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_3; x = 0, 0.4, 5, 10, 40\n";
    for (int i = 0; i < nq; ++i)
      {
	auto q = i * delq;
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {_Tp{0.0L}, _Tp{0.4L}, _Tp{5.0L}, _Tp{10.0L}, _Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_3(q, x);
	std::cout << '\n';
      }
    std::cout.flush();

    std::cout << "\n\n\n";
    std::cout << "theta_4; x = 0, 0.4, 5, 10, 40\n";
    for (int i = 0; i < nq; ++i)
      {
	auto q = i * delq;
	std::cout << ' ' << std::setw(w) << q;
	for (auto x : {_Tp{0.0L}, _Tp{0.4L}, _Tp{5.0L}, _Tp{10.0L}, _Tp{40.0L}})
	  std::cout << ' ' << std::setw(w) << std::__detail::__jacobi_theta_4(q, x);
	std::cout << '\n';
      }
    std::cout.flush();
  }


int
main()
{
  std::__detail::__jacobi_theta_1(0.205, 40.0);
  std::__detail::__jacobi_theta_1(0.210, 40.0);

  plot_jacobi_theta(1.0);

  test_jacobi_theta(1.0);
}
