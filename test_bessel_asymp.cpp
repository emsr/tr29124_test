/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_bessel_asymp test_bessel_asymp.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_bessel_asymp
*/

#include <iostream>
#include <sstream>
#include <ext/cmath>


  template <typename _Tp>
    void
    __cyl_bessel_jn_asymp_old(_Tp __nu, _Tp __x,
			      _Tp & __Jnu, _Tp & __Nnu)
    {
      const auto __2nu = 2 * __nu;
      const auto __x8 = 8 * __x;
      auto __k = 1;
      auto __k2m1 = 1;
      auto __P = _Tp(1);
      auto __Q = (__2nu - __k2m1) * (__2nu + __k2m1) / __x8;
      ++__k;
      const auto __eps = std::numeric_limits<_Tp>::epsilon();
      auto __t = _Tp(1);
      do
        {
          __k2m1 += 2;
          __t *= -(__2nu - __k2m1) * (__2nu + __k2m1) / (__k * __x8);
          auto __convP = std::abs(__t) < __eps * std::abs(__P);
          __P += __t;
          ++__k;

          __k2m1 += 2;
          __t *= (__2nu - __k2m1) * (__2nu + __k2m1) / (__k * __x8);
          auto __convQ = std::abs(__t) < __eps * std::abs(__Q);
          __Q += __t;
          ++__k;

          if (__convP && __convQ && __k > (__nu / _Tp(2)))
            break;
        }
      while (__k < 50 * __nu);

      auto __chi = __x - (__nu + _Tp(0.5L))
        	       * __gnu_cxx::__math_constants<_Tp>::__pi_half;
      auto __c = std::cos(__chi);
      auto __s = std::sin(__chi);

      auto __coef = std::sqrt(_Tp(2)
        	  / (__gnu_cxx::__math_constants<_Tp>::__pi * __x));
      __Jnu = __coef * (__c * __P - __s * __Q);
      __Nnu = __coef * (__s * __P + __c * __Q);

      return;
    }

int
main(int n_app_args, char ** app_arg)
{
  double x = 1000.0;
  double nu = 20;

  bool use_internal = false;
  if (n_app_args > 1)
  {
    int i_internal = 0;
    std::istringstream in(app_arg[1]);
    in >> i_internal;
    if (i_internal != 0)
      use_internal = true;
  }
  if (n_app_args > 2)
  {
    std::istringstream in(app_arg[2]);
    in >> nu;
  }
  if (n_app_args > 3)
  {
    std::istringstream in(app_arg[3]);
    in >> x;
  }

  std::cout << '\n';
  std::cout << "use internal function = " << std::boolalpha << use_internal << '\n';

  std::cout << '\n';
  std::cout << "nu = " << nu << '\n';

  do
  {
    __gnu_cxx::__cyl_bessel_t<double, double, double> Bess;
    try
    {
      Bess = std::__detail::__cyl_bessel_jn(nu, x);
    }
    catch (std::exception e)
    {
      std::cout << '\n' << "Couldn't run main Bessel function." << '\n';
    }

    __gnu_cxx::__cyl_bessel_t<double, double, double> BessAsym;
    if (use_internal)
      {
        double Jnua = 0.0, Nnua = 0.0, Jpnua = 0.0, Npnua = 0.0;
        __cyl_bessel_jn_asymp_old(nu, x, Jnua, Nnua);
        BessAsym.__nu_arg = nu;
        BessAsym.__x_arg = x;
        BessAsym.__J_value = Jnua;
        BessAsym.__J_deriv = Jpnua;
        BessAsym.__N_value = Nnua;
        BessAsym.__N_deriv = Npnua;
      }
    else
      BessAsym = std::__detail::__cyl_bessel_jn_asymp(nu, x);

    std::cout << '\n';
    std::cout << "x = " << x << '\n';
    std::cout << "Jnu = " << Bess.__J_value << '\n';
    std::cout << "Nnu = " << Bess.__N_value << '\n';
    std::cout << "Jnua = " << BessAsym.__J_value << '\n';
    std::cout << "Nnua = " << BessAsym.__N_value << '\n';
    std::cout << "Jnua - Jnu = " << BessAsym.__J_value - Bess.__J_value << '\n';
    std::cout << "Nnua - Nnu = " << BessAsym.__N_value - Bess.__N_value << '\n';
    if (!use_internal)
    {
      std::cout << "Jpnu = " << Bess.__J_deriv << '\n';
      std::cout << "Npnu = " << Bess.__N_deriv << '\n';
      std::cout << "Jpnua = " << BessAsym.__J_deriv << '\n';
      std::cout << "Npnua = " << BessAsym.__N_deriv << '\n';
      std::cout << "Jpnua - Jpnu = " << BessAsym.__J_deriv - Bess.__J_deriv << '\n';
      std::cout << "Npnua - Npnu = " << BessAsym.__N_deriv - Bess.__N_deriv << '\n';
    }

    x += 1000.0;
  }
  while (x < 100000.0);

  return 0;
}
