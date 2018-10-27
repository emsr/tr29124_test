/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_pow test_pow.cpp -lquadmath
./test_pow > test_pow.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_pow test_pow.cpp -lquadmath
./test_pow > test_pow.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * Integer powers.
 */
template<typename _Tp>
  _Tp
  pow(_Tp __x, int __n)
  {
    _Tp __val = _Tp{1};

    if (__n < 0)
      {
	__n = -__n;
	if (__x == _Tp{})
	  {
	    auto __u = _Tp{1} / __x;
	    // Correct sign of infinity.
	    __val = (__n % 2) ? __u : (__u * __u) ;
          }

	__x = _Tp{1} / __x;
      }

    // Repeated squaring method
    // Returns 0^0 = 1, so continuous in x.
    do
      {
	if ((__n & 1) == 1)
	  __val *= __x;
	__n >>= 1;
	__x *= __x;
      }
    while (__n > 0);

    return __val;
  }

int
main()
{
  for (int n = -5; n <= 5; ++n)
    {
      for (int i = 0; i <= 100; ++i)
	{
	  auto x = i * 0.1;
	  std::cout << ' ' << x
		    << ' ' << pow(x, n)
		    << ' ' << std::pow(x, n)
		    << '\n';
	}
    }
}
