// $HOME/bin/bin/g++ -o test_gamma test_gamma.cpp

// LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_gamma

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>

template<typename _Tp>
  _Tp
  a(_Tp g, unsigned int k)
  {
    
  }

template<typename _Tp>
  void
  lanczos(_Tp z, _Tp g)
  {
    // From pugh..
    int n = -2 - 0.3 * std::log(std::numeric_limits<_Tp>::epsilon);

    std::vector<_Tp> a;
    for (unsigned int k = 1; k < 100; ++k)
      {
	;
      }

    constexpr auto _S_log_sqrt_2pi = 9.189385332046727417803297364056176398620e-1L;
    auto logGammap1 = _S_log_sqrt_2pi + (z + 0.5L) * std::log(z + g + 0.5L)
	       - (z + g + 0.5L);
    auto fact = _Tp{1};
    auto sum = 0.5L * a(g, 0);
    for (unsigned int k = 1; k < 100; ++k)
      {
	fact *= (z - k - 1) / (z + k);
	auto term = fact * a(g, k);
      }
  }

template<typename _Tp>
  void
  spouge()
  {
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_2pi = _Tp{6.283185307179586476925286766559005768391L};
    auto a = _Tp{1};
    const auto __fact = _Tp{1} / std::sqrt(_S_2pi);
    while (_S_eps >= __fact * std::pow(_S_2pi, -a) / std::sqrt(a);
      a += _Tp{1};

    auto __term = _Tp{1};
    auto __sum = __term * std::sqrt(a - 1) * std::exp(a - 1);
    for (int __k = 2; __k < std::ceil(a); ++__k)
      {
	__term *= -_Tp{1} / _Tp(__k);
	__sum += __term * std::pow(a - __k, k - 0.5) * std::exp(a - __k);
      }
  }

int
main()
{
}
