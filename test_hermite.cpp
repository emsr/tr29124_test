// $HOME/bin_specfun/bin/g++ -std=gnu++1z -o test_hermite test_hermite.cpp -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_hermite > test_hermite.txt

#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <bits/float128.h>
#include "new_hermite.tcc"

int
main()
{
  std::cout.precision(std::numeric_limits<__float128>::digits10);
  std::cout.flags(std::ios::showpoint);
  auto width = 6 + std::cout.precision();

  std::cout << "\f\n\n  Table of integer sqrt\n";
  std::cout << "  =====================\n";
  for (int n = 0; n <= 100; ++n)
    std::cout << ' ' << std::setw(width) << std::sqrt(__float128(n)) << '\n';

  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout.flags(std::ios::showpoint);
  width = 6 + std::cout.precision();

  std::cout << "\f\n\n  Examine asymptotic transition region\n";
  std::cout << "  ====================================\n";
  for (int n = 4; n <= 200; ++n)
    {
      auto xt = std::sqrt(2.0 * n);
      auto del = 0.5 * xt / 200;
      std::cout << "  " << std::setw(width) << "n = " << std::setw(width) << n << '\n';
      std::cout << "  " << std::setw(width) << "x_t = " << std::setw(width) << xt << '\n';
      std::cout << "  " << std::setw(width) << "x";
      std::cout << "  " << std::setw(width) << "Hr_" << n << "(x)";
      std::cout << "  " << std::setw(width) << "Ha_" << n << "(x)";
      std::cout << '\n';
      for (int i = 0; i <= 201; ++i)
        {
          auto x = 0.75 * xt + i * del;
          auto h = __poly_hermite_recursion(n, x);
          auto ht = __poly_hermite_asymp(n, x);
	  if (std::abs(x - xt) < del)
	    std::cout << ">>";
	  else if (std::abs(x - 0.95 * xt) < del)
	    std::cout << "> ";
	  else if (std::abs(x - 1.05 * xt) < del)
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

  std::cout << "\f\n\n  Compare recursion and asymptotic\n";
  std::cout << "  ================================\n";
  for (int n = 0; n <= 50; ++n)
    {
      auto xt = std::sqrt(2.0 * n);
      std::cout << "  " << std::setw(width) << "n = " << std::setw(width) << n << '\n';;
      std::cout << "  " << std::setw(width) << "x_t = " << std::setw(width) << xt << '\n';;
      std::cout << "  " << std::setw(width) << "x";
      std::cout << "  " << std::setw(width) << "Hr_" << n << "(x)";
      std::cout << "  " << std::setw(width) << "Ha_" << n << "(x)";
      std::cout << '\n';
      auto del = 0.1;
      for (int i = 0; i <= 100; ++i)
        {
          auto x = i * del;
          auto h = __poly_hermite_recursion(n, x);
          auto ht = __poly_hermite_asymp(n, x);
	  if (std::abs(x - xt) < del)
	    std::cout << ">>";
	  else if (std::abs(x - 0.95 * xt) < del)
	    std::cout << "> ";
	  else if (std::abs(x - 1.05 * xt) < del)
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
  auto power = 1.0;
  for (int n = 0; n <= 50; ++n)
    {
      // sqrt(factorial(n) * 2**n sqrt(pi))
      auto norm = std::exp(0.5 * std::lgamma(n + 1)) * std::sqrt(power * std::sqrt(M_PI));
      power *= 2.0;
      std::cout << "  " << std::setw(width) << "n = " << n << '\n';
      std::cout << "  " << std::setw(width) << "norm = " << norm << '\n';
      std::cout << "  " << std::setw(width) << "x";
      std::cout << "  " << std::setw(width) << "H_" << n << "(x)";
      std::cout << "  " << std::setw(width) << "H~_" << n << "(x)";
      std::cout << "  " << std::setw(width) << "ratio";
      std::cout << '\n';
      for (int i = 0; i <= 100; ++i)
        {
          auto x = i * 0.1;
          auto h = __poly_hermite_recursion(n, x);
          auto ht = __poly_hermite_norm_recursion(n, x);
          std::cout << "  " << std::setw(width) << x;
          std::cout << "  " << std::setw(width) << h;
          std::cout << "  " << std::setw(width) << ht;
          std::cout << "  " << std::setw(width) << h / ht;
          std::cout << '\n';
        }
    }
  return 0;
}
