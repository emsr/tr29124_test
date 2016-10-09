/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -I. -o test_reperiodized_hyper test_reperiodized_hyper.cpp wrap_boost.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_reperiodized_hyper > test_reperiodized_hyper.txt

g++ -std=gnu++17 -g -DNO_SINH_COSH_PI -I. -o test_reperiodized_hyper test_reperiodized_hyper.cpp wrap_boost.cpp -lgsl -lgslcblas -lquadmath
./test_reperiodized_hyper > test_reperiodized_hyper.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>

#include <ext/cmath>
//#include <bits/float128.h>
#include <bits/sf_trig.tcc>

#include "wrap_boost.h"

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{

#ifdef NO_SINH_COSH_PI

  // Reperiodized sine function.

  /**
   * Return the reperiodized sine function @f$ \sinh_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see sinh_pi for more details.
   */
  inline float
  sinh_pif(float __x)
  { return std::__detail::__sinh_pi<float>(__x); }

  /**
   * Return the reperiodized sine function @f$ \sinh_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see sinh_pi for more details.
   */
  inline long double
  sinh_pil(long double __x)
  { return std::__detail::__sinh_pi<long double>(__x); }

  /**
   * Return the reperiodized hyperbolic sine function @f$ \sinh_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized hyperbolic sine function is defined by:
   * @f[
   * 	\sinh_\pi(x) = \sinh(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sinh_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__sinh_pi<__type>(__x);
    }

  // Reperiodized hyperbolic cosine function.

  /**
   * Return the reperiodized hyperbolic cosine function @f$ \cosh_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see cosh_pi for more details.
   */
  inline float
  cosh_pif(float __x)
  { return std::__detail::__cosh_pi<float>(__x); }

  /**
   * Return the reperiodized hyperbolic cosine function @f$ \cosh_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see cosh_pi for more details.
   */
  inline long double
  cosh_pil(long double __x)
  { return std::__detail::__cosh_pi<long double>(__x); }

  /**
   * Return the reperiodized hyperbolic cosine function @f$ \cosh_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized hyperbolic cosine function is defined by:
   * @f[
   * 	\cosh_\pi(x) = \cosh(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    cosh_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__cosh_pi<__type>(__x);
    }

#endif // NO_SINH_COSH_PI

  // Reperiodized hyperbolic tangent function.

  /**
   * Return the reperiodized hyperbolic tangent function @f$ \tanh_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see tanh_pi for more details.
   */
  inline float
  tanh_pif(float __x)
  { return std::__detail::__tanh_pi<float>(__x); }

  /**
   * Return the reperiodized hyperbolic tangent function @f$ \tanh_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see tanh_pi for more details.
   */
  inline long double
  tanh_pil(long double __x)
  { return std::__detail::__tanh_pi<long double>(__x); }

  /**
   * Return the reperiodized hyperbolic tangent function @f$ \tanh_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized hyperbolic tangent function is defined by:
   * @f[
   * 	\tanh_\pi(x) = \tanh(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    tanh_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__tanh_pi<__type>(__x);
    }


} // namespace __gnu_cxx

template<typename _Tp>
  void
  run_sin_cosh_pi()
  {
    constexpr _Tp _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 10 + std::cout.precision();

    std::cout << std::endl;
    std::cout << std::setw(width) << " x"
	      << std::setw(width) << " sinh_pi (GCC)"
	      << std::setw(width) << " delta sinh(pi x)"
	      << std::setw(width) << " cosh_pi (GCC)"
	      << std::setw(width) << " delta cosh(pi x)"
	      << std::setw(width) << " tanh_pi (GCC)"
	      << std::setw(width) << " delta tanh(pi x)"
	      << '\n';
    std::cout << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << '\n';
    for (int i = -1600; i <= +1600; ++i)
      {
	auto x = _Tp(0.1Q * i);
	auto sinh_pi_g = __gnu_cxx::sinh_pi(x);
	auto cosh_pi_g = __gnu_cxx::cosh_pi(x);
	auto tanh_pi_g = __gnu_cxx::tanh_pi(x);
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << sinh_pi_g
		  << ' ' << std::setw(width) << sinh_pi_g - std::sinh(_S_pi * x)
		  << ' ' << std::setw(width) << cosh_pi_g
		  << ' ' << std::setw(width) << cosh_pi_g - std::cosh(_S_pi * x)
		  << ' ' << std::setw(width) << tanh_pi_g
		  << ' ' << std::setw(width) << tanh_pi_g - std::tanh(_S_pi * x)
		  << '\n';
      }
    std::cout << std::endl;
  }


int
main()
{
  std::cout << "\ndouble\n=====\n\n";
  run_sin_cosh_pi<double>();

  std::cout << "\nlong double\n=====\n\n";
  run_sin_cosh_pi<long double>();
}
