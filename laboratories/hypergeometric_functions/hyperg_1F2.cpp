/**
 *
 */

#include <limits>
#include <cmath>
#include <stdexcept>

template<typename Tp>
  Tp
  hyperg_1f2(Tp a, Tp b, Tp c, Tp x, Tp& err)
  {
    constexpr auto s_max_iter = 200;
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    const auto s_huge = std::numeric_limits<Tp>::max() / Tp{1000};
    auto an = a;
    auto bn = b;
    auto cn = c;
    auto a0 = Tp{1};
    auto sum = Tp{1};
    auto n = 1;
    auto t = Tp{1};
    auto max = Tp{0};
    const auto stop = std::numeric_limits<Tp>::epsilon();

    do
      {
	if (an == Tp{0})
	  break;
	if (bn == Tp{0})
	  throw std::runtime_error("hyperg_1f2: series failed");
	if (cn == Tp{0})
	  throw std::runtime_error("hyperg_1f2: series failed");
	if (a0 > s_huge || n > s_max_iter)
	  throw std::runtime_error("hyperg_1f2: series failed");
	a0 *= (an * x) / (bn * cn * n);
	sum += a0;
	an += Tp{1};
	bn += Tp{1};
	cn += Tp{1};
	++n;
	auto z = std::abs(a0);
	if (z > max)
	  max = z;
	if (sum != Tp{0})
	  t = std::abs(a0 / sum);
	else
	  t = z;
      }
    while (t > stop);

    err = std::abs(s_eps * max / sum);

    return sum;
  }

#ifdef MAIN_1F2
int
main()
{
  double a = 3.0;
  double b = 2.0;
  double c = 1.5;
  double x = 2.5;
  double err = 0.0;
  hyperg_1f2(a, b, c, x, err);
}
#endif
