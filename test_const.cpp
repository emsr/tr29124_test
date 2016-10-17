// $HOME/bin/bin/g++ -std=gnu++17 -I. -o test_const test_const.cpp -lquadmath -lmpfr

// LD_LIBRARY_PATH=/$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_const > test_const.txt

#include <mpreal.h>
#include <bits/float128.h>
#include <ext/cmath>
#include <ext/math_const_mpreal.h>
#include <bits/numeric_limits_mpreal.h>

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <typeindex>
#include <string>

template<typename _Tp>
  void
  test_const(_Tp proto)
  {
    using cnst = __gnu_cxx::__math_constants<long double>;

    auto name{std::type_index{typeid(proto)}.name()};

    auto prec = __gnu_cxx::__digits10(proto);
    auto wd = 8 + prec;
    auto lprec = std::numeric_limits<long double>::digits10;
    auto lwd = 8 + lprec;

    std::cout << '\n';
    std::cout << "type                      : " << std::quoted(name) << '\n';
    std::cout << " 4 pi/3                   : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__4_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_4_pi(proto) << '\n';
    std::cout << " 4 pi/3                   : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__4_pi_div_3
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_4_pi_div_3(proto) << '\n';
    std::cout << " 2 pi                     : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__2_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_2_pi(proto) << '\n';
    std::cout << " pi                       : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_pi(proto) << '\n';
    std::cout << " pi/2                     : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__pi_half
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_pi_half(proto) << '\n';
    std::cout << " pi/3                     : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__pi_third
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_pi_third(proto) << '\n';
    std::cout << " pi/4                     : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__pi_quarter
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_pi_quarter(proto) << '\n';
    std::cout << " sqrt(pi)                 : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__root_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_root_pi(proto) << '\n';
    std::cout << " cbrt(pi)                 : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__cbrt_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_cbrt_pi(proto) << '\n';
    std::cout << " sqrt(pi/2)               : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__root_pi_div_2
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_root_pi_div_2(proto) << '\n';
    std::cout << " 1/ pi                    : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__one_div_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_one_div_pi(proto) << '\n';
    std::cout << " 2/ pi                    : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__two_div_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_two_div_pi(proto) << '\n';
    std::cout << " 2/ sqrt(pi)              : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__two_div_root_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_two_div_root_pi(proto) << '\n';
    std::cout << " pi^2/6                   : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__pi_sqr_div_6
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_pi_sqr_div_6(proto) << '\n';
    std::cout << " sqrt(2 pi)               : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__root_2_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_root_2_pi(proto) << '\n';
    std::cout << " ln(sqrt(2 pi))           : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__ln_root_2_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_ln_root_2_pi(proto) << '\n';
    std::cout << " rad2deg                  : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__deg
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_deg(proto) << '\n';
    std::cout << " deg2rad                  : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__rad
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_rad(proto) << '\n';
    std::cout << " Euler's number e         : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__e
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_e(proto) << '\n';
    std::cout << " 1/e                      : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__one_div_e
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_one_div_e(proto) << '\n';
    std::cout << " log_2(e)                 : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__log2_e
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_log2_e(proto) << '\n';
    std::cout << " log_10(e)                : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__log10_e
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_log10_e(proto) << '\n';
    std::cout << " ln(2)                    : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__ln_2
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_ln_2(proto) << '\n';
    std::cout << " ln(3)                    : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__ln_3
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_ln_3(proto) << '\n';
    std::cout << " ln(10)                   : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__ln_10
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_ln_10(proto) << '\n';
    std::cout << " log(pi)                  : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__ln_pi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_ln_pi(proto) << '\n';
    std::cout << " Euler-Mascheroni gamma_E : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__gamma_e
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_gamma_e(proto) << '\n';
    std::cout << " Golden Ratio phi         : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__phi
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_phi(proto) << '\n';
    std::cout << " Catalan's constant       : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__catalan
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_catalan(proto) << '\n';
    std::cout << " sqrt(2)                  : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__root_2
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_root_2(proto) << '\n';
    std::cout << " sqrt(3)                  : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__root_3
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_root_3(proto) << '\n';
    std::cout << " sqrt(3)                  : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__root_3_div_2
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_root_3_div_2(proto) << '\n';
    std::cout << " sqrt(5)                  : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__root_5
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_root_5(proto) << '\n';
    std::cout << " sqrt(7)                  : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__root_7
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_root_7(proto) << '\n';
    std::cout << " cbrt(3)                  : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__cbrt_3
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_cbrt_3(proto) << '\n';
    std::cout << "1/ sqrt(2)                : "
	      << std::setprecision(lprec) << std::setw(lwd) << cnst::__one_div_root_2
	      << std::setprecision(prec) << std::setw(wd) << __gnu_cxx::__const_one_div_root_2(proto) << '\n';
  }

int
main()
{
  test_const(1.0f);
  test_const(1.0);
  test_const(1.0l);
  test_const(1.0q);
  test_const(mpfr::mpreal(1, 256));
}
