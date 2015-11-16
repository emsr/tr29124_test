// $HOME/bin_specfun/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_airy test_airy.cpp gsl_wrap.cpp -lgsl -lgslcblas

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_airy

#include <iostream>
#include <iomanip>
#include <cmath>

#include "gsl_wrap.h"

double
chebyshev_eval(double a, double b, double * c, int m, double x)
{
  if ((x - a) * (x - b) > 0.0)
  {
    std::cerr << "x not in range in routine chebyshev_eval.";
    return 0.0;
  }

  double d = 0.0, dd = 0.0, y;
  double y2 = 2.0 * (y = (2.0 * x - a - b) / ( b - a ));
  for (int j = m - 1; j >= 1; --j)
  {
    double sv = d;
    d = y2 * d - dd + c[j];
    dd = sv;
  }
  return y * d - dd + 0.5 * c[0];
}

void
bessel_cheb(double x, double &gam1, double &gam2, double &gampl, double &gammi)
{
  double xx;
  static double c1[]
  {
    -1.142022680371168,
     0.0065165112670737,
     0.0003087090173086,
    -0.0000034706269649,
     0.0000000069437664,
     0.0000000000367765,
    -0.0000000000001356
  };
  static double c2[]
  {
    1.843740587300905,
   -0.0786528408447867,
    0.0012719271366546,
   -0.0000049717367042,
   -0.0000000331261198,
    0.0000000002423096,
   -0.0000000000001702,
   -0.00000000000000149
  };

  xx = 8.0 * x * x - 1.0;
  gam1 = chebyshev_eval( -1.0, +1.0, c1, 7, xx );
  gam2 = chebyshev_eval( -1.0, +1.0, c2, 8, xx );
  gampl = gam2 - x * gam1;
  gammi = gam2 + x * gam1;
}

  template<typename _Tp>
    void
    __gamma_temme(_Tp __mu,
		  _Tp & __gam1, _Tp & __gam2, _Tp & __gampl, _Tp & __gammi)
    {
      constexpr _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr _Tp _S_gamma_E = std::__detail::__numeric_constants<_Tp>::__gamma_e();
      __gampl = _Tp{1} / std::tgamma(_Tp{1} + __mu);
      __gammi = _Tp{1} / std::tgamma(_Tp{1} - __mu);

      if (std::abs(__mu) < _S_eps)
	__gam1 = -_S_gamma_E;
      else
	__gam1 = (__gammi - __gampl) / (_Tp{2} * __mu);

      __gam2 = (__gammi + __gampl) / _Tp{2};

      return;
    }

int
main()
{
  std::cout.precision(16);
  std::cout.flags(std::ios::showpoint);

  for (auto i = -200; i <= +200; ++i)
    {
      auto x = i * 0.05;
      auto Ai = __gnu_cxx::airy_ai(x);
      auto Bi = __gnu_cxx::airy_bi(x);
      auto Ai_gsl = wrap_gsl_sf_airy_ai(x);
      auto Bi_gsl = wrap_gsl_sf_airy_bi(x);
      //double Ai_tr1, Bi_tr1, Aip_tr1, Bip_tr1;
      //std::tr1::__detail::__airy(x, Ai_tr1, Bi_tr1, Aip_tr1, Bip_tr1);
      std::cout << ' ' << std::setw(8) << x
                << ' ' << std::setw(20) << Ai
                << ' ' << std::setw(20) << Bi
                << ' ' << std::setw(20) << Ai_gsl
                << ' ' << std::setw(20) << Bi_gsl
                << ' ' << std::setw(20) << std::abs(Ai - Ai_gsl)
                << ' ' << std::setw(20) << std::abs(Bi - Bi_gsl) << '\n';
    }

  for (auto i = -10; i <= +10; ++i)
    {
      auto x = i * 0.05;
      double gam1, gam2, gampl, gammi;
      bessel_cheb(x, gam1, gam2, gampl, gammi);
      std::cout << '\n';
      std::cout << ' ' << std::setw(8) << x
                << ' ' << std::setw(20) << gam1
                << ' ' << std::setw(20) << gam2
                << ' ' << std::setw(20) << gampl
                << ' ' << std::setw(20) << gammi << '\n';
      __gamma_temme(x, gam1, gam2, gampl, gammi);
      std::cout << ' ' << std::setw(8) << x
                << ' ' << std::setw(20) << gam1
                << ' ' << std::setw(20) << gam2
                << ' ' << std::setw(20) << gampl
                << ' ' << std::setw(20) << gammi << '\n';
    }
}
