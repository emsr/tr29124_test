// $HOME/bin/bin/g++ -std=gnu++14 -o airy_toy_old airy_toy_old.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./airy_toy_old > airy_toy_old.txt

// g++ -std=gnu++14 -DNO_LOGBQ -o airy_toy_old airy_toy_old.cpp -L$HOME/bin/lib64 -lquadmath

// ./airy_toy_old > airy_toy_old.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>

#include "float128.h"

template<typename _Tp>
  void
  run_toy()
  {
    constexpr _Tp _F_k[9]
    {
      1.0 / 6.0,       //  1 / 3!==6
      4.0 / 720.0,     //  1*4 / 6!==720
      28.0 / 362880.0, //  1*4*7 / 9!==362880
      28.0 / 47900160.0,  //  (1*4*7*10 / 10) / (12! / 10)==47900160
      280.0 * 13.0 / 130767.0e7,  //  1*4*7*10*13 / (15!==1307674368000)
      280.0 * 208.0 / 640237.0e10,  //  1*4*7*10*13*16 / 18!==6402373705728000
      280.0 * 208.0 * 19.0 / 5.1090942e19,  //  1*4*7*10*13*16*19 / 21!==51090942171709440000
      280.0 * 208.0 * 418.0 / 6.2044840173e23,  //  1*4*7*10*13*16*19*22 / 24!==620448401733239439360000
      280.0 * 208.0 * 418.0 * 25.0 / 1.0888869e28  //  1*4*7*10*13*16*19*22*25 / 27!==10888869450418352160768000000
    };

    _Tp _Fp_k[9];
    for (int __k = 0; __k < 9; ++__k)
      _Fp_k[__k] = 3 * (__k + 1) * _F_k[__k];

    constexpr _Tp _G_k[9]
     {
      2.0 / 24.0,  //  2 / 4!
      10.0 / 5040.0,  //  2*5 / 7!
      80.0 / 36288.0e2,  //  2*5*8 / 10!==3628800
      880.0 / 622702.0e4,  //  2*5*8*11 / 13!==6227020800
      880.0 * 14.0 / 209228.0e8,  //  2*5*8*11*14 / 16!==20922789888000
      880.0 * 238.0 / 121645.0e12,  //  2*5*8*11*14*17 / 19!==121645100408832000
//            4760.0 was 4780.0!
      880.0 * 4760.0 / 1.1240007278e21,  //  2*5*8*11*14*17*20 / 22!==1124000727777607680000
      880.0 * 4760.0 * 23.0 / 1.551121e25,  //  2*5*8*11*14*17*20*23 / 25!==15511210043330985984000000
      880.0 * 4760.0 * 598.0 / 3.048883446e29  //  2*5*8*11*14*17*20*23*26 / 28!==304888344611713860501504000000
    };

    _Tp _Gp_k[9];
    for (int __k = 0; __k < 9; ++__k)
      _Gp_k[__k] = (3 * (__k + 1) + 1) * _G_k[__k];

    constexpr int _N_cd = 200;
    _Tp _S_cn[_N_cd], _S_dn[_N_cd];
    _S_cn[0] = _Tp{15} / _Tp{216};
    _S_dn[0] = -(_Tp{7} / _Tp{5}) * _S_cn[0];
    for (int __n = 2; __n <= _N_cd; ++__n)
      {
	_S_cn[__n - 1] = _S_cn[__n - 2]
		       * (6 * __n - 5) * (6 * __n - 3) * (6 * __n - 1)
		       / (216 * __n * (2 * __n - 1));
	_S_dn[__n - 1] = -_S_cn[__n - 1] * (6 * __n + 1) / (6 * __n - 1);
      }

    std::cout.precision(std::numeric_limits<_Tp>::max_digits10);
    std::cout << "\nc\n";
    for (const auto& c : _S_cn)
      std::cout << c << '\n';
    std::cout << "\nd\n";
    for (const auto& d : _S_dn)
      std::cout << d << '\n';
    std::cout << '\n';
    std::cout << '\n';
    std::cout << "\nF\n";
    std::cout << 1 << '\n';
    for (const auto& c : _F_k)
      std::cout << c << '\n';
    std::cout << "\nFp\n";
    std::cout << 1 << '\n';
    for (const auto& c : _Fp_k)
      std::cout << c << '\n';
    std::cout << "\nG\n";
    std::cout << 1 << '\n';
    for (const auto& c : _G_k)
      std::cout << c << '\n';
    std::cout << "\nGp\n";
    std::cout << 1 << '\n';
    for (const auto& c : _Gp_k)
      std::cout << c << '\n';
  }

int
main()
{
  run_toy<long double>();
  run_toy<__float128>();
}
