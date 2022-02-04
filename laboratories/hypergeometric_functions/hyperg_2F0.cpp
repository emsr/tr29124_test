/**
 *
 */

#include <limits>
#include <cmath>
#include <stdexcept>

template<typename Tp>
  Tp
  hyperg_2f0(Tp a, Tp b, Tp x, Tp& err)
  {
    constexpr auto s_max_iter = 200;
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    const auto s_huge = std::numeric_limits<Tp>::max() / Tp{1000};
    auto an = a;
    auto bn = b;
    auto a0 = Tp{1};
    auto sum = Tp{1};
    auto n = 1;
    auto t = Tp{1};
    auto max = Tp{0};
    auto conv = std::numeric_limits<Tp>::max() / Tp{2};
    auto conv1 = conv;
    const auto stop = std::numeric_limits<Tp>::epsilon();

    do
      {
	if (an == Tp{0})
	  break;
	if (bn == Tp{0})
	  break;
	if (a0 > s_huge || n > s_max_iter)
	  throw std::runtime_error("hyperg_2f0: series failed");
	a0 *= (an * bn * x) / n;
	an += Tp{1};
	bn += Tp{1};
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

    t = std::abs(s_eps * max / sum);
    max = std::abs(conv / sum);
    if (max > t)
      t = max;
    err = t;

    return sum;
  }

#ifdef MAIN_2F0
int
main()
{
  double a = -3.0;
  double b = 2.0;
  double x = 2.5;
  double err = 0.0;
  hyperg_2f0(a, b, x, err);
}
#endif

