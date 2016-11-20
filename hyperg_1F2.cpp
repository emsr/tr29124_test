/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o hyperg_1F2 hyperg_1F2.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./hyperg_1F2 > hyperg_1F2.txt
*/

#include <limits>
#include <cmath>

template<typename _Tp>
  _Tp
  __hyperg_1f2(_Tp a, _Tp b, _Tp c, _Tp x, _Tp& err)
  {
    constexpr auto _S_max_iter = 200;
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_huge = std::numeric_limits<_Tp>::max() / _Tp{1000};
    auto an = a;
    auto bn = b;
    auto cn = c;
    auto a0 = _Tp{1};
    auto sum = _Tp{1};
    auto n = 1;
    auto t = _Tp{1};
    auto max = _Tp{0};

    do
      {
	if (an == _Tp{0})
	  break;
	if (bn == _Tp{0})
	  throw std::runtime_error("__hyperg_1f2: series failed");
	if (cn == _Tp{0})
	  throw std::runtime_error("__hyperg_1f2: series failed");
	if (a0 > _S_huge || n > _S_max_iter)
	  throw std::runtime_error("__hyperg_1f2: series failed");
	a0 *= (an * x) / (bn * cn * n);
	sum += a0;
	an += _Tp{1};
	bn += _Tp{1};
	cn += _Tp{1};
	++n;
	auto z = std::abs(a0);
	if (z > max)
	  max = z;
	if (sum != _Tp{0})
	  t = std::abs(a0 / sum);
	else
	  t = z;
      }
    while (t > stop);

    err = std::abs(_S_eps * max / sum);

    return sum;
  }

int
main()
{
  double a = 3.0;
  double b = 2.0;
  double c = 1.5;
  double x = 2.5;
  double err = 0.0;
  __hyperg_1f2(a, b, c, x, err);
}
