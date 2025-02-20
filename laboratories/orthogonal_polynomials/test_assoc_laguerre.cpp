/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>

#include <emsr/special_functions.h>

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
			<< ' ' << std::setw(12) << emsr::assoc_laguerre(n, alpha, x)
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
			<< ' ' << std::setw(12) << emsr::assoc_laguerre(n, alpha, x)
			<< '\n';
	    }
	}
    }
}
