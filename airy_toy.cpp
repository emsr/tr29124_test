// $HOME/bin/bin/g++ -std=gnu++14 -o airy_toy airy_toy.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./airy_toy

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

    _S_cn.push_back(_Tp{15} / _Tp{216});
    _S_dn.push_back(-(_Tp{7} / _Tp{5}) * _S_cn.back());
    for (int __n = 2; __n <= 200; ++__n)
      {
	_S_cn.push_back(_S_cn.back()
		       * (6 * __n - 5) * (6 * __n - 3) * (6 * __n - 1)
		       / (216 * __n * (2 * __n - 1)));
	_S_dn.push_back(-_S_cn.back() * (6 * __n + 1) / (6 * __n - 1));
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
    for (int __k = 0; __k < 9; ++__k)
      {
	_F_k.push_back();
	_Fp_k.push_back(3 * (__k + 1) * _F_k.back());
	_G_k.push_back();
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
