/*
$HOME/bin_tr29124/bin/g++ -Wall -Wextra -o test_conf_hyperg_limit test_conf_hyperg_limit.cpp wrap_gsl.cpp -lgsl -lgslcblas 2> err.txt
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_conf_hyperg_limit > test_conf_hyperg_limit.txt

g++ -std=c++14 -DNO_CBRT -DNO_LOGBQ -Wall -Wextra -o test_conf_hyperg_limit test_conf_hyperg_limit.cpp wrap_gsl.cpp -lgsl -lgslcblas 2> err.txt
./test_conf_hyperg_limit > test_conf_hyperg_limit.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

template<typename _Tp>
  _Tp
  conf_hyperg_limit_sum(_Tp __c, _Tp __z)
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
  conf_hyperg_limit(_Tp __c, _Tp __z)
  {
    return conf_hyperg_limit_sum(__c, __z);
  }

int
main()
{
  constexpr auto prec = std::numeric_limits<double>::digits10;
  std::cout.precision(prec);
  double c = 0.2;
  for (int i = -200; i < +200; ++i)
  {
    double z = 0.1 * i;
    std::cout << ' ' << std::setw(6) << z << ' ' << conf_hyperg_limit(c, z) << '\n';
  }
}
