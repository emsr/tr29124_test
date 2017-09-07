/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I.. -o test_gauss_hermite test_gauss_hermite.cpp -lquadmath
*/

#include <cmath>

#include "gauss_hermite_integrate.h"

template<typename _Tp>
  struct my_f
  {
    _Tp a;
    _Tp b;
    _Tp c;
    _Tp d;
    _Tp e;

    _Tp
    operator()(_Tp x)
    { return (((a * x + b) * x + c ) * x + d ) * x + e; }
  };

template<typename _Tp>
  struct my_fc
  {
    _Tp a;
    _Tp b;
    _Tp c;

    _Tp
    operator()(_Tp x)
    { return a * std::cos(b * x + c); }
  };

template<typename _Tp>
  void
  test_hermite_quad(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    std::cout << '\n';

    my_f<_Tp> params{ 0.5, 0.0, 0.0, -1.0, 1.3 };
    my_fc<_Tp> cparams{ 1.0, 1.0, 0.0 };

    //pfoo2 = __gnu_ext::gauss_hermite_prob_integrate<_Tp>(params, 3);
    //std::cout << "integral = "  << std::setw(width)<< pfoo2 << '\n';
    //std::cout << "delta = " << std::setw(width) << pfoo2 - _Tp{7.01855916896680140676414279747L} << '\n';
    //pfooc = __gnu_ext::gauss_hermite_prob_integrate<_Tp>(cparams, 29);
    //std::cout << "integral = " << std::setw(width) << pfooc << '\n';
    //std::cout << "delta = "  << std::setw(width)<< pfooc - _Tp{1.52034690106628080561194014675L} << '\n';

    auto xfoo2 = __gnu_ext::gauss_hermite_integrate<_Tp>(params, 3);
    std::cout << "integral = " << std::setw(width) << xfoo2 << '\n';
    std::cout << "delta = " << std::setw(width) << xfoo2 - _Tp{2.96886020026673934572443053460L} << '\n';
    auto xfooc = __gnu_ext::gauss_hermite_integrate<_Tp>(cparams, 29);
    std::cout << "integral = " << std::setw(width) << xfooc << '\n';
    std::cout << "delta = " << std::setw(width) << xfooc - _Tp{1.38038844704314297477341524673L} << '\n';
  }

int
main()
{
  test_hermite_quad<double>();
  test_hermite_quad<long double>();
}
