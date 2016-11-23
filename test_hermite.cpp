/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_hermite test_hermite.cpp -lquadmath
./test_hermite > test_hermite.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_hermite test_hermite.cpp -lquadmath
./test_hermite > test_hermite.txt
*/

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <bits/float128_io.h>
#include "new_hermite.tcc"

//#include "quadrature/integration.h"

template<typename _Tp>
  void
  test_hermite(_Tp proto)
  {
    const auto _S_pi = __gnu_cxx::__const_pi(proto);

    auto fname = [](std::string_view front, int n, std::string_view back)
		 {
		   std::ostringstream out;
		   out << front << n << back;
		   return out.str();
		 };

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    std::cout << "\f\n\n  Table of integer sqrt\n";
    std::cout << "  =====================\n";
    for (int n = 0; n <= 100; ++n)
      std::cout << "  " << std::setw(4) << n
		<< "  " << std::setw(width) << std::sqrt(_Tp(n)) << '\n';

    std::cout << "\f\n\n  Table of factorial sqrt\n";
    std::cout << "  =======================\n";
    auto fact = _Tp{1};
    for (int n = 0; n <= 100; ++n)
      {
	std::cout << "  " << std::setw(4) << n
		  << "  " << std::setw(width) << fact << '\n';
	fact *= std::sqrt(_Tp(n + 1));
      }

    std::cout << "\f\n\n  Examine asymptotic transition region\n";
    std::cout << "  ====================================\n";
    for (int n = 4; n <= 200; ++n)
      {
	auto xt = std::sqrt(_Tp{2} * n);
	auto del = _Tp{0.5Q} * xt / _Tp{200};
	std::cout << "  " << std::setw(6) << "n   = " << std::setw(4) << n << '\n';
	std::cout << "  " << std::setw(6) << "x_t = " << std::setw(width) << xt << '\n';
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("Hr_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Ha_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta";
	std::cout << "  " << std::setw(width) << "delta-ratio";
	std::cout << "  " << std::setw(width) << fname("Ha2_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta2";
	std::cout << "  " << std::setw(width) << "delta2-ratio";
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	for (int i = 0; i <= 201; ++i)
          {
            auto x = _Tp{0.75Q} * xt + i * del;
            auto h = __poly_hermite_recursion(n, x);
            auto ht = __poly_hermite_asymp(n, x);
            auto h2 = __poly_hermite_asymp2(n, x);
	    if (std::abs(x - xt) < del)
	      std::cout << ">>";
	    else if (std::abs(x - _Tp{0.95Q} * xt) < del)
	      std::cout << "> ";
	    else if (std::abs(x - _Tp{1.05Q} * xt) < del)
	      std::cout << "> ";
            else
              std::cout << "  ";
	    std::cout << std::setw(width) << x
		      << "  " << std::setw(width) << h
		      << "  " << std::setw(width) << ht
		      << "  " << std::setw(width) << ht - h
		      << "  " << std::setw(width) << (ht - h) / h
		      << "  " << std::setw(width) << h2
		      << "  " << std::setw(width) << h2 - h
		      << "  " << std::setw(width) << (h2 - h) / h
		      << '\n';
          }
      }

    std::cout << "\f\n\n  Examine asymptotic transition region\n";
    std::cout << "  ====================================\n";
    for (int n : {1000, 2000, 5000, 10000})
      {
	auto xt = std::sqrt(_Tp{2} * n);
	auto del = _Tp{0.5Q} * xt / _Tp{200};
	std::cout << "  " << std::setw(6) << "n   = " << std::setw(4) << n << '\n';
	std::cout << "  " << std::setw(6) << "x_t = " << std::setw(width) << xt << '\n';
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("Hr_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Ha_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta";
	std::cout << "  " << std::setw(width) << "delta-ratio";
	std::cout << "  " << std::setw(width) << fname("Ha2_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta2";
	std::cout << "  " << std::setw(width) << "delta2-ratio";
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	for (int i = 1; i <= 201; ++i)
          {
            auto x = _Tp{0.75Q} * xt + i * del;
            auto h = __poly_hermite_recursion(n, x);
            auto ht = __poly_hermite_asymp(n, x);
            auto h2 = __poly_hermite_asymp2(n, x);
	    if (std::abs(x - xt) < del)
	      std::cout << ">>";
	    else if (std::abs(x - _Tp{0.95Q} * xt) < del)
	      std::cout << "> ";
	    else if (std::abs(x - _Tp{1.05Q} * xt) < del)
	      std::cout << "> ";
            else
              std::cout << "  ";
	    std::cout << std::setw(width) << x
		      << "  " << std::setw(width) << h
		      << "  " << std::setw(width) << ht
		      << "  " << std::setw(width) << ht - h
		      << "  " << std::setw(width) << (ht - h) / h
		      << "  " << std::setw(width) << h2
		      << "  " << std::setw(width) << h2 - h
		      << "  " << std::setw(width) << (h2 - h) / h
		      << '\n';
          }
      }

    std::cout << "\f\n\n  Compare recursion and asymptotic\n";
    std::cout << "  ================================\n";
    for (int n = 0; n <= 50; ++n)
      {
	auto xt = std::sqrt(_Tp{2} * n);
	std::cout << "  " << std::setw(6) << "n   = " << std::setw(4) << n << '\n';;
	std::cout << "  " << std::setw(6) << "x_t = " << std::setw(width) << xt << '\n';;
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("Hr_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Ha_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta";
	std::cout << "  " << std::setw(width) << "delta-ratio";
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	auto del = _Tp{0.1Q};
	for (int i = 0; i <= 100; ++i)
          {
            auto x = i * del;
            auto h = __poly_hermite_recursion(n, x);
            auto ht = __poly_hermite_asymp(n, x);
	    if (std::abs(x - xt) < del)
	      std::cout << ">>";
	    else if (std::abs(x - _Tp{0.95Q} * xt) < del)
	      std::cout << "> ";
	    else if (std::abs(x - _Tp{1.05Q} * xt) < del)
	      std::cout << "> ";
            else
              std::cout << "  ";
            std::cout << std::setw(width) << x
		      << "  " << std::setw(width) << h
		      << "  " << std::setw(width) << ht
		      << "  " << std::setw(width) << ht - h
		      << "  " << std::setw(width) << (ht - h) / h
		      << '\n';
          }
      }

    std::cout << "\f\n\n  Compare normalizations\n";
    std::cout << "  ======================\n";
    auto power = _Tp{1};
    for (int n = 0; n <= 50; ++n)
      {
	// sqrt(factorial(n) * 2**n sqrt(pi))
	auto norm = std::exp(_Tp{0.5Q} * std::lgamma(n + 1)) * std::sqrt(power * std::sqrt(_S_pi));
	power *= _Tp{2};
	std::cout << "  " << std::setw(width) << "n = " << n << '\n';
	std::cout << "  " << std::setw(width) << "norm = " << norm << '\n';
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("H_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Hn_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("He_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Hen_", n, "(x)");
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	for (int i = 0; i <= 100; ++i)
          {
            auto x = i * _Tp{0.1Q};
            auto h = __poly_hermite_recursion(n, x);
            auto hn = __poly_hermite_norm_recursion(n, x);
            auto he = __poly_prob_hermite_recursion(n, x);
            auto hen = __poly_prob_hermite_norm_recursion(n, x);
            std::cout << "  " << std::setw(width) << x;
            std::cout << "  " << std::setw(width) << h;
            std::cout << "  " << std::setw(width) << hn;
            std::cout << "  " << std::setw(width) << he;
            std::cout << "  " << std::setw(width) << hen;
            std::cout << '\n';
          }
      }
    return;
  }

template<typename _Tp>
  void
  test_hermite_norm(_Tp proto)
  {
    return;
  }


int
main()
{
  test_hermite(1.0Q);
  //test_hermite(1.0F);
}
