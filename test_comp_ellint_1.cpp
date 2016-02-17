// $HOME/bin/bin/g++ -std=c++1z -o test_comp_ellint_1 test_comp_ellint_1.cpp

// LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_comp_ellint_1 > test_comp_ellint_1.txt

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <experimental/utility>

auto IMAX = 0;

//  Use AGM to do an ab initio calculation of K(k).
template<typename _Tp>
  _Tp
  __ellint_K(_Tp __k)
  {
    constexpr auto _S_max_iter = 100;
    auto __am = _Tp{0.5} * (_Tp{1} + __k);
    auto __gm = std::sqrt(__k);
    for (int __i = 0; __i < _S_max_iter; ++__i)
      {
IMAX = __i;
	__gm = std::sqrt(__gm * std::exchange(__am, _Tp{0.5} * (__am + __gm)));
	if (std::abs(__am - __gm) < std::numeric_limits<_Tp>::epsilon())
	  break;
      }
    return __gm;
  }

//  Use AGM to do an ab initio calculation of the elliptic nome.
template<typename _Tp>
  _Tp
  __ellint_nome(_Tp __k)
  {
    constexpr auto _S_pi = _Tp{3.1415926535897932384626433832795029L};
    auto __kp = std::sqrt((_Tp{1} - __k) * (_Tp{1} + __k));
    auto __K = __ellint_K(__k);
    auto __Kp = __ellint_K(__kp);
    return std::exp(-_S_pi * __Kp / __K);
  }

template<typename _Tp>
  void
  test_K()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);

    auto width = 6 + std::cout.precision();
    for (int i = 0; i <= 1000; ++i)
      {
	auto k = _Tp(i * 0.001L);
	std::cout << ' ' << IMAX
		  << ' ' << std::setw(width) << k
		  << ' ' << std::setw(width) << __ellint_K(k)
		  << ' ' << std::setw(width) << __ellint_nome(k) << '\n';
      }
  }

int
main()
{
  std::cout << "\n\nfloat\n";
  test_K<float>();

  std::cout << "\n\ndouble\n";
  test_K<double>();

  std::cout << "\n\nlong double\n";
  test_K<long double>();
}
