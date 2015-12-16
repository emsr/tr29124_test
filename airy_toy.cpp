// $HOME/bin/bin/g++ -std=gnu++14 -o airy_toy airy_toy.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./airy_toy > airy_toy.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>

#include "float128.h"

template<typename _Tp>
  void
  run_toy()
  {
    std::vector<_Tp> _S_cn, _S_dn;

    _S_cn.push_back(_Tp{1});
    _S_dn.push_back(-_Tp{1});
    auto __denom = _Tp{1};
    for (int __s = 1; __s <= 200; ++__s)
      {
	//  This also works actually.
	//_S_cn.push_back(_S_cn.back()
	//	       * (6 * __s - 5) * (6 * __s - 3) * (6 * __s - 1)
	//	       / (216 * __s * (2 * __s - 1)));
	__denom *= 216 * __s;
	auto __numer = _Tp{1};
	for (int __r = 0; __r < 2 * __s; ++__r)
	  __numer *= (2 * __s + 2 * __r + 1);
	_S_cn.push_back(__numer / __denom);
	_S_dn.push_back(-_S_cn.back() * (6 * __s + 1) / (6 * __s - 1));
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
    std::vector<_Tp> _F_k, _Fp_k, _G_k, _Gp_k;
    _F_k.push_back(_Tp{1});
    _Fp_k.push_back(_Tp{1});
    _G_k.push_back(_Tp{1});
    _Gp_k.push_back(_Tp{1});
    auto __numer = _Tp{1};
    auto __denom = _Tp{1};
    for (int __k = 1; __k < 9; ++__k)
      {
	__numer *= (3 * __k - 2) * (3 * __k - 1) * (3 * __k - 0);
	__denom *= __k * (__k + 1) * (__k + 2);
	_F_k.push_back(_Tp{1} / __numer);
	_Fp_k.push_back(3 * (__k + 1) * _F_k.back());
	_G_k.push_back(_Fp_k.back());
	_Gp_k.push_back((3 * (__k + 1) + 1) * _G_k.back());
      }
    std::cout << "\nF\n";
    for (const auto& c : _F_k)
      std::cout << c << '\n';
    std::cout << "\nFp\n";
    for (const auto& c : _Fp_k)
      std::cout << c << '\n';
    std::cout << "\nG\n";
    for (const auto& c : _G_k)
      std::cout << c << '\n';
    std::cout << "\nGp\n";
    for (const auto& c : _Gp_k)
      std::cout << c << '\n';
  }

int
main()
{
  run_toy<long double>();
  run_toy<__float128>();
}
