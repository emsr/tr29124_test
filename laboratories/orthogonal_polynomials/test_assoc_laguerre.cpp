/*
$HOME/bin/bin/g++ -std=c++2a -g -Wall -Wextra -Wno-psabi -I../../include -I../../cxx_fp_utils/include -I../../polynomial/include -I../../quadrature/include -I../../cxx_summation/include -o test_assoc_laguerre test_assoc_laguerre.cpp
*/

#include <cmath>
#include <iostream>
#include <iomanip>

int
main()
{
  for (unsigned int n : {0, 1, 2, 5})
    {
      for (double alpha : {0.0, 1.0/3.0, 0.5, 2.0/3.0, 1.0})
	{
	  std::cout << "\n\n n = " << std::setw(2) << n
		    << "  alpha = " << std::setw(12) << alpha << '\n';
	  for (int i = -100; i <= 200; ++i)
	    {
	      double x = 0.1 * i;
	      std::cout << ' ' << std::setw(12) << x
			<< ' ' << std::setw(12) << std::assoc_laguerre(n, alpha, x)
			<< '\n';
	    }
	}
    }

  for (unsigned int n : {0, 1, 2, 5})
    {
      for (int alpha : {-1, -2, -3})
	{

	  std::cout << "\n\n n = " << std::setw(2) << n
		    << "  alpha = " << std::setw(12) << alpha << '\n';
	  for (int i = -100; i <= 200; ++i)
	    {
	      double x = 0.1 * i;
	      std::cout << ' ' << std::setw(12) << x
			<< ' ' << std::setw(12) << std::assoc_laguerre(n, alpha, x)
			<< '\n';
	    }
	}
    }
}
