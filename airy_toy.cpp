// $HOME/bin/bin/g++ -std=gnu++14 -o airy_toy airy_toy.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./airy_toy > airy_toy.txt

// g++ -std=gnu++14 -DNO_LOGBQ -o airy_toy airy_toy.cpp -lquadmath

// ./airy_toy > airy_toy.txt

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
    for (int __s = 1; __s <= 200; ++__s)
      {
        // Turn this into a recursion:
	// for (int r = 0; r < 2 * s; ++r)
	//   numer *= (2 * s + 2 * r + 1);
	_S_cn.push_back(_S_cn.back()
		      * (6 * __s - 5) * (6 * __s - 3) * (6 * __s - 1)
		      / (216 * __s * (2 * __s - 1)));
	_S_dn.push_back(-_S_cn.back() * (6 * __s + 1) / (6 * __s - 1));
      }

    std::cout.precision(std::numeric_limits<_Tp>::max_digits10);
    std::cout << std::showpoint;

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
    auto __F_numer = _Tp{1};
    auto __F_denom = _Tp{1};
    auto __G_numer = _Tp{1};
    auto __G_denom = _Tp{1};
    for (unsigned long long __k = 1ULL; __k <= 15ULL; ++__k)
      {
	std::cout << '\n' << ' ' << std::setw(40) << __F_numer << '/' << std::setw(40) << __F_denom;
	__F_denom *= (3ULL * __k - 2ULL) * (3ULL * __k - 1ULL) * (3ULL * __k - 0ULL);
	__F_numer *= 1ULL + 3ULL * (__k - 1ULL);
	_F_k.push_back(__F_numer / __F_denom);
	_Fp_k.push_back(3 * __k * _F_k.back());

	std::cout << '\t' << ' ' << std::setw(40) << __G_numer << '/' << std::setw(40) << __G_denom;
	__G_denom *= (3ULL * __k - 1ULL) * (3ULL * __k) * (3ULL * __k + 1ULL);
	__G_numer *= 2ULL + 3ULL * (__k - 1ULL);
	_G_k.push_back(__G_numer / __G_denom);
	_Gp_k.push_back((3 * __k + 1) * _G_k.back());
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
  std::cout << "\nfloat\n-----\n";
  run_toy<float>();

  std::cout << "\ndouble\n------\n";
  run_toy<double>();

  std::cout << "\nlong double\n-----------\n";
  run_toy<long double>();

  std::cout << "\n__float128\n----------\n";
  run_toy<__float128>();
}
