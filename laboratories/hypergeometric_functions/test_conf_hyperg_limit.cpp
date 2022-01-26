/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

  template<typename _Tp>
    _Tp
    conf_hyperg_limit_sum(_Tp c, _Tp z)
    {
      constexpr int s_max_iter = 10000;
      _Tp term{1};
      _Tp sum = term;
      for (int i = 0; i < s_max_iter; ++i)
	{
	  term *=  z / ((c + i) * (i + 1));
	  sum += term;
	  if (std::abs(term) < std::numeric_limits<_Tp>::epsilon())
	    break;
	}
      return sum;
    }

  template<typename _Tp>
    _Tp
    conf_hyperg_limit(_Tp c, _Tp z)
    {
      return conf_hyperg_limit_sum(c, z);
    }

template<typename _Tp>
  void
  test_conf_hyperg_limit(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto c = _Tp{0.2Q};
    const auto del = _Tp{1} / _Tp{10};
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
