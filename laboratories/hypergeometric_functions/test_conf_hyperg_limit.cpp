/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

#include <emsr/numeric_limits.h>
#include <emsr/sf_hyperg.h>

  template<typename Tp>
    Tp
    conf_hyperg_limit_sum(Tp c, Tp z)
    {
      constexpr int s_max_iter = 10000;
      Tp term{1};
      Tp sum = term;
      for (int i = 0; i < s_max_iter; ++i)
	{
	  term *=  z / ((c + i) * (i + 1));
	  sum += term;
	  if (std::abs(term) < std::numeric_limits<Tp>::epsilon())
	    break;
	}
      return sum;
    }

  template<typename Tp>
    Tp
    conf_hyperg_limit(Tp c, Tp z)
    {
      return conf_hyperg_limit_sum(c, z);
    }

template<typename Tp>
  void
  test_conf_hyperg_limit(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto c = Tp{0.2Q};
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i < +200; ++i)
    {
      auto z = del * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << emsr::conf_hyperg_lim(c, z)
		<< '\n';
    }
  }

int
main()
{
  test_conf_hyperg_limit(1.0);
}
