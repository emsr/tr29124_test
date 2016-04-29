
// $HOME/bin_specfun/bin/g++ -g -std=c++14 -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o carlson_std carlson_std.cpp

#include <ext/cmath>
#include <iostream>
#include <complex>
#include <limits>

int
main()
{
  using namespace std::literals::complex_literals;
  std::cout.precision(std::numeric_limits<double>::digits10);

  std::cout << "R_F(1, 2, 0)         = " << __gnu_cxx::ellint_rf(1.0, 2.0, 0.0) << std::endl;  //  1.3110287771461
  std::cout << "R_F(i, -i, 0)        = " << __gnu_cxx::ellint_rf(1.0i, -1.0i, 0.0) << std::endl;  //  1.8540746773014;
  std::cout << "R_F(i - 1, i, 0)     = " << __gnu_cxx::ellint_rf(-1.0 + 1.0i, 1.0i, 0.0) << std::endl;  //  0.79612586584234 - i 1.2138566698365
  std::cout << "R_F(0.5, 1, 0)       = " << __gnu_cxx::ellint_rf(0.5, 1.0, 0.0) << std::endl;  //  1.8540746773014
  std::cout << "R_F(2, 3, 4)         = " << __gnu_cxx::ellint_rf(2.0, 3.0, 4.0) << std::endl;  //  0.58408284167715
  std::cout << "R_F(i, -i, 2)        = " << __gnu_cxx::ellint_rf(1.0i, -1.0i, 2.0) << std::endl;  //  1.0441445654064
  std::cout << "R_F(i - 1, 1 - i, i) = " << __gnu_cxx::ellint_rf(-1.0 + 1.0i, 1.0 - 1.0i, 1.0i) << std::endl;  //  0.93912050218619 - i 0.53296252018635

  std::cout << std::endl;

  std::cout << "R_C(0, 1/4)          = " << __gnu_cxx::ellint_rc(0.0, 0.25) << std::endl;  //  pi = 3.1415926535898
  std::cout << "R_C(9/4, 2)          = " << __gnu_cxx::ellint_rc(2.25, 2.0) << std::endl;  //  ln 2 = 0.69314718055995
  std::cout << "R_C(0, i)            = " << __gnu_cxx::ellint_rc(0.0, 1.0i) << std::endl;  //  (1 - i) 1.1107207345396
  std::cout << "R_C(-i, i)           = " << __gnu_cxx::ellint_rc(-1.0i, 1.0i) << std::endl;  //  1.2260849569072 - i 0.34471136988768
  std::cout << "R_C(1/4, -2)         = " << __gnu_cxx::ellint_rc(0.25, -2.0) << std::endl;  //  ln 2/3 = 0.23104906018665
  std::cout << "R_C(i, -1)           = " << __gnu_cxx::ellint_rc(1.0i, -1.0) << std::endl;  //  0.77778596920447 + i 0.19832484993429

  std::cout << std::endl;

  std::cout << "R_J(0, 1, 2, 3)      = " << __gnu_cxx::ellint_rj(0.0, 1.0, 2.0, 3.0) << std::endl;  //  0.77688623778582
  std::cout << "R_J(2, 3, 4, 5)      = " << __gnu_cxx::ellint_rj(2.0, 3.0, 4.0, 5.0) << std::endl;  //  0.14297579667157
  std::cout << "R_J(2, 3, 4, -1 + i) = " << __gnu_cxx::ellint_rj(2.0, 3.0, 4.0, -1.0 + 1.0i) << std::endl;  //  0.13613945827771 - i 0.38207561624427
  std::cout << "R_J(i, -i, 0, 2)     = " << __gnu_cxx::ellint_rj(1.0i, -1.0i, 0.0, 2.0) << std::endl;  //  1.6490011662711
  std::cout << "R_J(-1 + i, -1 - i, 1, 2) = " << __gnu_cxx::ellint_rj(-1.0 + 1.0i, -1.0 - 1.0i, 1.0, 2.0) << std::endl;  //  0.94148358841220
  std::cout << "R_J(i, -i, 0, 1 - i)      = " << __gnu_cxx::ellint_rj(1.0i, -1.0i, 0.0, 1.0 - 1.0i) << std::endl;  //  1.8260115229009 + i 1.2290661908643
  std::cout << "R_J(-1 + i, -1 - i, 1, -3 + i)  = " << __gnu_cxx::ellint_rj(-1.0 + 1.0i, -1.0 - 1.0i, 1.0, -3.0 + 1.0i) << std::endl;  //  -0.61127970812028 - i 1.0684038390007
  std::cout << "R_J(-1 + i, -2 - i, -i, -1 + i) = " << __gnu_cxx::ellint_rj(-1.0 + 1.0i, -2.0 - 1.0i, 0.0 - 1.0i, -1.0 + 1.0i) << std::endl;  //  1.8249027393704 - i 1.2218475784827

  std::cout << std::endl;

  std::cout << "R_D(0, 2, 1)            = " << __gnu_cxx::ellint_rd(0.0, 2.0, 1.0) << std::endl;  //  1.7972103521034
  std::cout << "R_D(2, 3, 4)            = " << __gnu_cxx::ellint_rd(2.0, 3.0, 4.0) << std::endl;  //  0.16510527294261
  std::cout << "R_D(i, -i, 2)           = " << __gnu_cxx::ellint_rd(1.0i, -1.0i, 2.0) << std::endl;  //  0.65933854154220
  std::cout << "R_D(0, i, -i)           = " << __gnu_cxx::ellint_rd(0.0, 1.0i, -1.0i) << std::endl;  //  1.2708196271910 + i 2.7811120159521
  std::cout << "R_D(0, i - 1, i)        = " << __gnu_cxx::ellint_rd(0.0, -1.0 + 1.0i, 1.0i) << std::endl;  //  -1.8577235439239 - i 0.96193450888839
  std::cout << "R_D(-2 - i, -i, -1 + i) = " << __gnu_cxx::ellint_rd(-2.0 - 1.0i, -1.0i, -1.0 + 1.0i) << std::endl;  //  1.8249027393704 - i 1.2218475784827

  std::cout << std::endl;

  std::cout << "R_G(0, 16, 16)    = 2E(0) = pi = " << __gnu_cxx::ellint_rg(0.0, 16.0, 16.0) << std::endl;  //  3.1415926535898
  std::cout << "R_G(2, 3, 4)      = " << __gnu_cxx::ellint_rg(2, 3, 4) << std::endl;  //  1.7255030280692
  std::cout << "R_G(0, i, -i)     = " << __gnu_cxx::ellint_rg(0.0, 1.0i, -1.0i) << std::endl;  //  0.42360654239699
  std::cout << "R_G(i - 1, i, 0)  = " << __gnu_cxx::ellint_rg(-1.0 + 1.0i, 1.0i, 0.0) << std::endl;  //  0.44660591677018 + i 0.70768352357515
  std::cout << "R_G(-i, i - 1, i) = " << __gnu_cxx::ellint_rg(-1.0i, -1.0 + 1.0i, 1.0i) << std::endl;  //  0.36023392184473 + i 0.40348623401722
  std::cout << "R_G(0, 0.0796, 4) = E(0.99) = " << __gnu_cxx::ellint_rg(0, 0.0796, 4) << std::endl;  //  1.0284758090288

  std::cout << std::endl;

  std::cout << "R_F<int>(1, 2, 0) = " << __gnu_cxx::ellint_rf<int>(1.0, 2.0, 0.0) << std::endl;  //  1.3110287771461

  std::cout << std::endl;

  double __pi_2 = __gnu_cxx::__math_constants<double>::__pi_half;
  std::cout << "K(pi/2)  = " << std::__detail::__comp_ellint_1(__pi_2 + 0.0i) << std::endl;  //  1.5887715763658593607082818553065 - i 1.3986463677643598308560440635658
  std::cout << "K(-1)    = " << std::__detail::__comp_ellint_1(-1.0 + 0.0i) << std::endl;  //  1.31102877714605990523242
  std::cout << "K(2)     = " << std::__detail::__comp_ellint_1(2.0 + 0.0i) << std::endl;  //  1.31102877714605990523242 - i 1.31102877714605990523242

  return 0;
}
