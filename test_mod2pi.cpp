/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_mod2pi test_mod2pi.cpp -lquadmath
./test_mod2pi > test_mod2pi.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_mod2pi test_mod2pi.cpp -lquadmath
./test_mod2pi > test_mod2pi.txt
*/

#include <limits>
#include <iostream>
#include <ext/cmath>

template<typename _Tp>
  _Tp
  __mod2pi_cheap(_Tp __x)
  {
    const auto _S_2pi = __gnu_cxx::__const_2_pi(__x);
    return __x - _S_2pi * std::floor(__x / _S_2pi);
  }

template<typename _Tp>
  _Tp
  __mod2pi_cephes_wtf(_Tp __x)
  {
    const auto _S_2pi = __gnu_cxx::__const_2_pi(__x);
    const auto __n = std::floor(__x / _S_2pi);
    a = z - ldexp(n, 2);  /* 4n */
    a -= ldexp( n, 1);    /* 2n */
    a -= ldexp( n, -2 );  /* n/4 */
    a -= ldexp( n, -5 );  /* n/32 */
    a -= ldexp( n, -9 );  /* n/512 */
    a += ldexp( n, -15 ); /* add n/32768 */
    a -= ldexp( n, -17 ); /* n/131072 */
    a -= ldexp( n, -18 );
    a -= ldexp( n, -20 );
    a -= ldexp( n, -22 );
    a -= ldexp( n, -24 );
    a -= ldexp( n, -28 );
    a -= ldexp( n, -32 );
    a -= ldexp( n, -37 );
    a -= ldexp( n, -39 );
    a -= ldexp( n, -40 );
    a -= ldexp( n, -42 );
    a -= ldexp( n, -46 );
    a -= ldexp( n, -47 );

  }

template<typename _Tp>
  void
  test_mod2pi(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__max_digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;
  }

int
main()
{
  test_mod2pi<float>();
  test_mod2pi<double>();
  test_mod2pi<long double>();
  test_mod2pi<__float128>();
  //test_mod2pi<mpfr::mpreal>((1,  256));
}
