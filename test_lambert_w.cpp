/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_lambert_w test_lambert_w.cpp -lquadmath
./test_lambert_w > test_lambert_w.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_lambert_w test_lambert_w.cpp -lquadmath
./test_lambert_w > test_lambert_w.txt
*/

#include <ext/cmath>

  /**
   * This is the third-order Halley's root finding algorithm for Lambert W.
   */
  template<typename _Tp>
    _Tp
    __lambert_w_iter(_Tp __z)
    {
      const auto _S_max_iter = 1000u;

      for (auto __k = 0u; __k < _S_max_iter; ++__k)
	{
          auto __expwk = std::exp(__wk);
	  auto __fact = __wk * __expwk - __z
          __wkp1 -= __fact
		 / ((__wk + 1) * __expwk - (__wk + 2) * __fact / (2 * __wk + 2));
	  __wk = __wkp1;
	}
    }


template<typename _Tp>
  void
  test_lambert_w(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

  }

int
main()
{
  test_lambert_w(1.0);
}
