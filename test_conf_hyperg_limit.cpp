/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_conf_hyperg_limit test_conf_hyperg_limit.cpp -lquadmath -L. -lwgsl -lburkhardt
./test_conf_hyperg_limit > test_conf_hyperg_limit.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_conf_hyperg_limit test_conf_hyperg_limit.cpp -lquadmath -L. -lwgsl -lburkhardt
./test_conf_hyperg_limit > test_conf_hyperg_limit.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

  template<typename _Tp>
    _Tp
    __conf_hyperg_limit_sum(_Tp __c, _Tp __z)
    {
      constexpr int _S_max_iter = 10000;
      _Tp __term{1};
      _Tp __sum = __term;
      for (int __i = 0; __i < _S_max_iter; ++__i)
	{
	  __term *=  __z / ((__c + __i) * (__i + 1));
	  __sum += __term;
	  if (std::abs(__term) < std::numeric_limits<_Tp>::epsilon())
	    break;
	}
      return __sum;
    }

  template<typename _Tp>
    _Tp
    __conf_hyperg_limit(_Tp __c, _Tp __z)
    {
      return conf_hyperg_limit_sum(__c, __z);
    }

template<typename _Tp>
  void
  test_conf_hyperg_limit(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto c = _Tp{0.2Q};
    const auto del = _Tp{1} / _Tp{10};
    for (int i = -200; i < +200; ++i)
    {
      auto z = del * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << __gnu_cxx::conf_hyperg_lim(c, z)
		<< '\n';
    }
  }

int
main()
{
  test_conf_hyperg_limit(1.0);
}
