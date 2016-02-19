// $HOME/bin_specfun/bin/g++ -std=gnu++14 -o airy_toy airy_toy.cpp -L$HOME/bin/lib64 -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./airy_toy > airy_toy.txt

// g++ -std=gnu++14 -DNO_LOGBQ -I. -o airy_toy airy_toy.cpp -lquadmath

// ./airy_toy > airy_toy.txt

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <ext/math_const.h>
#include <bits/float128.h>
#include <bits/numeric_limits.h>

template<typename _Tp>
  void
  run_toy()
  {
    using __cmplx = std::complex<_Tp>;

    constexpr auto _S_eps = __gnu_cxx::__epsilon(_Tp{});
    /*constexpr*/ auto _S_log10min = __gnu_cxx::__log10_min(_Tp{});
    constexpr auto _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
    constexpr auto _S_sqrt_pi = __gnu_cxx::__math_constants<_Tp>::__root_pi;
    constexpr auto _S_Ai0{3.550280538878172392600631860041831763980e-1};
    constexpr auto _S_Aip0{2.588194037928067984051835601892039634793e-1};
    //constexpr auto _S_Bi0{6.149266274460007351509223690936135535960e-1};
    //constexpr auto _S_Bip0{8.868776642045783582807775119976424596506e-1};
    constexpr auto _S_1d6 = _Tp{1} / _Tp{6};
    constexpr auto _S_i = __cmplx{_Tp{0}, _Tp{1}};
    constexpr auto _S_big = _Tp{5.0L}; // was 3.5
    constexpr int _N_FG = 40;

    std::vector<_Tp> _S_cn, _S_dn;
    _S_cn.push_back(_Tp{1});
    _S_dn.push_back(-_Tp{1});
    for (int __s = 1; __s <= 200; ++__s)
      {
        // Turn this into a recursion:
	// for (int r = 0; r < 2 * s; ++r)
	//   numer *= (2 * s + 2 * r + 1);
	//auto __a = _S_cn.back()
	//	 * (6 * __s - 5) * (6 * __s - 3) * (6 * __s - 1)
	//	      / (216 * __s * (2 * __s - 1));
	auto __a = _S_cn.back()
		 * (_Tp(__s - 1) / _Tp{2} + _Tp{5} / _Tp(72 * __s));
	auto __b = -__a * _Tp(__s + _S_1d6) / _Tp(__s - _S_1d6);
	if (std::isnan(__a) || std::isinf(__a)
	 || std::isnan(__b) || std::isinf(__b))
	  break;
	_S_cn.push_back(__a);
	_S_dn.push_back(__b);
      }

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;

    std::cout << "\nc[" << _S_cn.size() << "]\n";
    for (const auto& c : _S_cn)
      std::cout << c << '\n';
    std::cout << "\nd[" << _S_dn.size() << "]\n";
    for (const auto& d : _S_dn)
      std::cout << d << '\n';

    std::cout << "\nc_even\n";
    for (int __s = _S_cn.size() - 1; __s >= 0; --__s)
      if (__s % 2 == 0)
	std::cout << _S_cn[__s] << '\n';
    std::cout << "\nc_odd\n";
    for (int __s = _S_cn.size() - 1; __s >= 0; --__s)
      if (__s % 2 == 1)
	std::cout << _S_cn[__s] << '\n';
    std::cout << "\nd_even\n";
    for (int __s = _S_dn.size() - 1; __s >= 0; --__s)
      if (__s % 2 == 0)
	std::cout << _S_dn[__s] << '\n';
    std::cout << "\nd_odd\n";
    for (int __s = _S_dn.size() - 1; __s >= 0; --__s)
      if (__s % 2 == 1)
	std::cout << _S_dn[__s] << '\n';

    std::cout << '\n';
    std::cout << '\n';
    std::vector<_Tp> _Fai, _Faip, _Gai, _Gaip, _Hai, _Haip;
    _Fai.push_back(_Tp{1});
    _Faip.push_back(_Tp{1});
    _Gai.push_back(_Tp{1});
    _Gaip.push_back(_Tp{1});
    _Hai.push_back(_Tp{1});
    _Haip.push_back(_Tp{1});
    auto __Fai_numer = _Tp{1};
    auto __Fai_denom = _Tp{1};
    auto __Gai_numer = _Tp{1};
    auto __Gai_denom = _Tp{1};
    auto __Hai_numer = _Tp{1};
    auto __Hai_denom = _Tp{1};
    const auto _S_min = std::numeric_limits<_Tp>::min();
    const auto __k_max = 40ULL;
    for (unsigned long long __k = 1ULL; __k <= __k_max; ++__k)
      {
	std::cout << '\n' << ' ' << std::setw(40) << __Fai_numer << '/' << std::setw(40) << __Fai_denom;
	//__Fai_denom *= (3ULL * __k - 2ULL) * (3ULL * __k - 1ULL) * (3ULL * __k);
	//__Fai_numer *= 1ULL + 3ULL * (__k - 1ULL);
	__Fai_denom *= (3ULL * __k - 1ULL) * (3ULL * __k);
	__Fai_numer *= 1ULL;
	if (__Fai_numer / __Fai_denom < _Tp{10} * _S_min)
	  break;
	_Fai.push_back(__Fai_numer / __Fai_denom);
	_Faip.push_back(3ULL * __k * _Fai.back());

	std::cout << '\t' << ' ' << std::setw(40) << __Gai_numer << '/' << std::setw(40) << __Gai_denom;
	//__Gai_denom *= (3ULL * __k - 1ULL) * (3ULL * __k) * (3ULL * __k + 1ULL);
	//__Gai_numer *= 2ULL + 3ULL * (__k - 1ULL);
	__Gai_denom *= (3ULL * __k) * (3ULL * __k + 1ULL);
	__Gai_numer *= 1ULL;
	_Gai.push_back(__Gai_numer / __Gai_denom);
	_Gaip.push_back((3ULL * __k + 1ULL) * _Gai.back());

	std::cout << '\t' << ' ' << std::setw(40) << __Hai_numer << '/' << std::setw(40) << __Hai_denom;
	//__Hai_denom *= (3ULL * __k) * (3ULL * __k + 1ULL) * (3ULL * __k + 2ULL);
	//__Hai_numer *= 3ULL + 3ULL * (__k - 1ULL);
	__Hai_denom *= (3ULL * __k + 1ULL) * (3ULL * __k + 2ULL);
	__Hai_numer *= 1;
	_Hai.push_back(__Hai_numer / __Hai_denom);
	_Haip.push_back((3ULL * __k + 2ULL) * _Hai.back());
      }
    std::cout << '\n';

    std::cout << "\nF[" << _Fai.size() << "]\n";
    for (const auto& c : _Fai)
      std::cout << c << '\n';
    std::cout << "\nFp[" << _Faip.size() << "]\n";
    for (const auto& c : _Faip)
      std::cout << c << '\n';
    std::cout << "\nG[" << _Gai.size() << "]\n";
    for (const auto& c : _Gai)
      std::cout << c << '\n';
    std::cout << "\nGp[" << _Gaip.size() << "]\n";
    for (const auto& c : _Gaip)
      std::cout << c << '\n';
    std::cout << "\nH[" << _Hai.size() << "]\n";
    for (const auto& c : _Hai)
      std::cout << c << '\n';
    std::cout << "\nHp[" << _Haip.size() << "]\n";
    for (const auto& c : _Haip)
      std::cout << c << '\n';
  }

int
main()
{
  std::cout << "\nfloat\n=====\n";
  run_toy<float>();

  std::cout << "\ndouble\n======\n";
  run_toy<double>();

  std::cout << "\nlong double\n===========\n";
  run_toy<long double>();

  std::cout << "\n__float128\n==========\n";
  run_toy<__float128>();
}
