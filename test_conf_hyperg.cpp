/*
$HOME/bin_tr29124/bin/g++ -g -Wall -Wextra -Wno-psabi -I. -o test_conf_hyperg test_conf_hyperg.cpp wrap_gsl.cpp $HOME/tr29124_test/gslextras/Fresnel/fresnel.c $HOME/tr29124_test/gslextras/Jacobi/jacobi-0.9.2/src/jacobi.c $HOME/tr29124_test/gslextras/Hermite/gsl_sf_hermite.c -lgsl -lgslcblas -lquadmath
./test_conf_hyperg > test_conf_hyperg.txt

g++ -std=c++14 -DNO_CBRT -DNO_LOGBQ -g -Wall -Wextra -I. -o test_conf_hyperg test_conf_hyperg.cpp wrap_gsl.cpp -lgsl -lgslcblas
./test_conf_hyperg > test_conf_hyperg.txt
*/

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

template<typename _Tp>
  void
  test_conf_hyperg(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto c = _Tp{0.2Q};
    for (int i = -200; i < +200; ++i)
    {
      auto z = _Tp{0.1Q} * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << std::conf_hyperg(a, c, z)
		<< '\n';
    }
  }

int
main()
{
}
