// $HOME/bin_specfun/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_bessel_asymp test_bessel_asymp.cpp

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_bessel_asymp


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
  double x = 100.0;
  double nu = 20;

  bool use_internal = false;
  bool do_spot = false;
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
    do_spot = true;
    std::istringstream in(app_arg[3]);
    in >> x;
  }

  std::cout << '\n';
  std::cout << "use internal function = " << std::boolalpha << use_internal << '\n';

  std::cout << '\n';
  std::cout << "nu = " << nu << '\n';

  do
  {
    double Jnu = 0.0, Nnu = 0.0, Jpnu = 0.0, Npnu = 0.0;
    try
    {
      std::__detail::__bessel_jn(nu, x, Jnu, Nnu, Jpnu, Npnu);
    }
    catch (std::exception e)
    {
      std::cout << '\n' << "Couldn't run main Bessel function." << '\n';
    }
    double Jnua = 0.0, Nnua = 0.0, Jpnua = 0.0, Npnua = 0.0;
    if (use_internal)
      __cyl_bessel_jn_asymp_old(nu, x, Jnua, Nnua);
    else
      std::__detail::__cyl_bessel_jn_asymp(nu, x, Jnua, Nnua, Jpnua, Npnua);
    std::cout << '\n';
    std::cout << "x = " << x << '\n';
    std::cout << "Jnu = " << Jnu << '\n';
    std::cout << "Nnu = " << Nnu << '\n';
    std::cout << "Jnua = " << Jnua << '\n';
    std::cout << "Nnua = " << Nnua << '\n';
    std::cout << "Jnua - Jnu = " << Jnua - Jnu << '\n';
    std::cout << "Nnua - Nnu = " << Nnua - Nnu << '\n';
    if (!use_internal)
    {
      std::cout << "Jpnu = " << Jpnu << '\n';
      std::cout << "Npnu = " << Npnu << '\n';
      std::cout << "Jpnua = " << Jpnua << '\n';
      std::cout << "Npnua = " << Npnua << '\n';
      std::cout << "Jpnua - Jpnu = " << Jpnua - Jpnu << '\n';
      std::cout << "Npnua - Npnu = " << Npnua - Npnu << '\n';
    }

    x += 1000.0;
  }
  while (x < 500.0 * nu);

  return 0;
}
