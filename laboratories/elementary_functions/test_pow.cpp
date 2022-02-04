/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>

/**
 * Integer powers.
 */
template<typename Tp>
  Tp
  pow(Tp x, int n)
  {
    Tp val = Tp{1};

    if (n < 0)
      {
	n = -n;
	if (x == Tp{})
	  {
	    auto u = Tp{1} / x;
	    // Correct sign of infinity.
	    val = (n % 2) ? u : (u * u) ;
          }

	x = Tp{1} / x;
      }

    // Repeated squaring method
    // Returns 0^0 = 1, so continuous in x.
    do
      {
	if ((n & 1) == 1)
	  val *= x;
	n >>= 1;
	x *= x;
      }
    while (n > 0);

    return val;
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
