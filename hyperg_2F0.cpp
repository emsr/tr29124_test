/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o hyperg_2F0 hyperg_2F0.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./hyperg_2F0 > hyperg_2F0.txt
*/

#include <limits>
#include <cmath>

template<typename _Tp>
  _Tp
  __hyperg_2f0(_Tp a, _Tp b, _Tp x, _Tp& err)
  {
    constexpr auto _S_max_iter = 200;
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_huge = std::numeric_limits<_Tp>::max() / _Tp{1000};
    auto an = a;
    auto bn = b;
    auto a0 = _Tp{1};
    auto sum = _Tp{1};
    auto n = 1;
    auto t = _Tp{1};
    auto max = _Tp{0};
    auto conv = std::numeric_limits<_Tp>::max() / _Tp{2};
    auto conv1 = conv;

    do
      {
	if (an == _Tp{0})
	  break;
	if (bn == _Tp{0})
	  break;
	if (a0 > _S_huge || n > _S_max_iter)
	  throw std::runtime_error("__hyperg_2f0: series failed");
	a0 *= (an * bn * x) / n;
	an += _Tp{1};
	bn += _Tp{1};
	++n;
	auto z = std::abs(a0);
	if (z > max)
	  max = z;
	if (z >= conv)
	  {
	    if (z < max && z > conv1)
	      break;
	  }
	conv1 = conv;
	conv = z;
	sum += a0;
	if (sum != 0)
	  t = std::abs(a0 / sum);
	else
	  t = z;
      }
    while (t > stop);

    t = std::abs(_S_eps * max / sum);
    max = std::abs(conv / sum);
    if (max > t)
      t = max;
    err = t;

    return sum;
  }

int
main()
{
  double a = -3.0;
  double b = 2.0;
  double x = 2.5;
  double err = 0.0;
  __hyperg_2f0(a, b, x, err);
}

