// $HOME/bin_specfun/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_bessel_asymp test_bessel_asymp.cpp

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_bessel_asymp


#include <iostream>
#include <sstream>
#include <cmath>


using std::__detail::__numeric_constants;

    template <typename _Tp>
    void
    __cyl_bessel_jn_asymp(const _Tp __nu, const _Tp __x,
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
        	       * __numeric_constants<_Tp>::__pi_2();
      auto __c = std::cos(__chi);
      auto __s = std::sin(__chi);

      auto __coef = std::sqrt(_Tp(2)
        	  / (__numeric_constants<_Tp>::__pi() * __x));
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

  std::cout << std::endl;
  std::cout << "use internal function = " << std::boolalpha << use_internal << std::endl;

  std::cout << std::endl;
  std::cout << "nu = " << nu << std::endl;

  do
  {
    double Jnu = 0.0, Nnu = 0.0, Jpnu = 0.0, Npnu = 0.0;
    try
    {
      std::__detail::__bessel_jn(nu, x, Jnu, Nnu, Jpnu, Npnu);
    }
    catch (std::exception e)
    {
      std::cout << std::endl << "Couldn't run main Bessel function." << std::endl;
    }
    double Jnua, Nnua;
    if (use_internal)
      __cyl_bessel_jn_asymp(nu, x, Jnua, Nnua);
    else
      std::__detail::__cyl_bessel_jn_asymp(nu, x, Jnua, Nnua);
    std::cout << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "Jnu = " << Jnu << std::endl;
    std::cout << "Nnu = " << Nnu << std::endl;
    std::cout << "Jnua = " << Jnua << std::endl;
    std::cout << "Nnua = " << Nnua << std::endl;
    std::cout << "Jnua - Jnu = " << Jnua - Jnu << std::endl;
    std::cout << "Nnua - Nnu = " << Nnua - Nnu << std::endl;

    x += 1000.0;
  }
  while (x < 500.0 * nu);

  return 0;
}
