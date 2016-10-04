/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -I. -o test_reperiodized_trig test_reperiodized_trig.cpp wrap_boost.cpp
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_reperiodized_trig > test_reperiodized_trig.txt

g++ -std=c++14 -o test_reperiodized_trig test_reperiodized_trig.cpp wrap_boost.cpp -lgsl -lgslcblas
./test_reperiodized_trig > test_reperiodized_trig.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include "wrap_boost.h"
#include <bits/sf_trig.tcc>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{

#ifdef NO_SIN_COS_PI

  // Reperiodized sine function.

  /**
   * Return the reperiodized sine function @f$ \sin_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see sin_pi for more details.
   */
  inline float
  sin_pif(float __x)
  { return std::__detail::__sin_pi<float>(__x); }

  /**
   * Return the reperiodized sine function @f$ \sin_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see sin_pi for more details.
   */
  inline long double
  sin_pil(long double __x)
  { return std::__detail::__sin_pi<long double>(__x); }

  /**
   * Return the reperiodized sine function @f$ \sin_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized sine function is defined by:
   * @f[
   * 	\sin_\pi(x) = \sin(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    sin_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__sin_pi<__type>(__x);
    }

  // Reperiodized cosine function.

  /**
   * Return the reperiodized cosine function @f$ \cos_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see cos_pi for more details.
   */
  inline float
  cos_pif(float __x)
  { return std::__detail::__cos_pi<float>(__x); }

  /**
   * Return the reperiodized cosine function @f$ \cos_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see cos_pi for more details.
   */
  inline long double
  cos_pil(long double __x)
  { return std::__detail::__cos_pi<long double>(__x); }

  /**
   * Return the reperiodized cosine function @f$ \cos_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized cosine function is defined by:
   * @f[
   * 	\cos_\pi(x) = \sin(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    cos_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__cos_pi<__type>(__x);
    }

#endif // NO_SIN_COS_PI

  // Reperiodized tangent function.

  /**
   * Return the reperiodized tangent function @f$ \tan_\pi(x) @f$
   * for @c float argument @f$ x @f$.
   *
   * @see tan_pi for more details.
   */
  inline float
  tan_pif(float __x)
  { return std::__detail::__tan_pi<float>(__x); }

  /**
   * Return the reperiodized tangent function @f$ \tan_\pi(x) @f$
   * for <tt>long double</tt> argument @f$ x @f$.
   *
   * @see tan_pi for more details.
   */
  inline long double
  tan_pil(long double __x)
  { return std::__detail::__tan_pi<long double>(__x); }

  /**
   * Return the reperiodized tangent function @f$ \tan_\pi(x) @f$
   * for real argument @f$ x @f$.
   *
   * The reperiodized tangent function is defined by:
   * @f[
   * 	\tan_\pi(x) = \sin(\pi x)
   * @f]
   *
   * @tparam _Tp The floating-point type of the argument @c __x.
   * @param __x The argument
   */
  template<typename _Tp>
    inline typename __gnu_cxx::__promote<_Tp>::__type
    tan_pi(_Tp __x)
    {
      typedef typename __gnu_cxx::__promote<_Tp>::__type __type;
      return std::__detail::__tan_pi<__type>(__x);
    }


} // namespace __gnu_cxx

template<typename _Tp>
  void
  run_sin_cos_pi()
  {
    constexpr _Tp _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << std::endl;
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sin_pi (GCC)"
	      << std::setw(width) << "sin_pi (Boost)"
	      << std::setw(width) << "delta sin_pi"
	      << std::setw(width) << "delta sin(pi x)"
	      << std::setw(width) << "cos_pi (GCC)"
	      << std::setw(width) << "cos_pi (Boost)"
	      << std::setw(width) << "delta cos_pi"
	      << std::setw(width) << "delta cos(pi x)"
	      << std::setw(width) << "tan_pi (GCC)"
	      << std::setw(width) << "delta tan(pi x)"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    for (int i = -1600; i <= +1600; ++i)
      {
	auto x = _Tp(0.1Q * i);
	auto sin_pi_g = __gnu_cxx::sin_pi(x);
	auto sin_pi_b = beast::sin_pi(x);
	auto cos_pi_g = __gnu_cxx::cos_pi(x);
	auto cos_pi_b = beast::cos_pi(x);
	auto tan_pi_g = __gnu_cxx::tan_pi(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << sin_pi_g
		  << std::setw(width) << sin_pi_b
		  << std::setw(width) << sin_pi_g - sin_pi_b
		  << std::setw(width) << sin_pi_g - std::sin(_S_pi * x)
		  << std::setw(width) << cos_pi_g
		  << std::setw(width) << cos_pi_b
		  << std::setw(width) << cos_pi_g - cos_pi_b
		  << std::setw(width) << cos_pi_g - std::cos(_S_pi * x)
		  << std::setw(width) << tan_pi_g
		  << std::setw(width) << tan_pi_g - std::tan(_S_pi * x)
		  << '\n';
      }
    std::cout << std::endl;

    std::cout << std::endl;
    std::cout << std::setw(width) << "x"
	      << std::setw(width) << "sin_pi (GCC)"
	      << std::setw(width) << "sin_pi (Boost)"
	      << std::setw(width) << "delta sin_pi"
	      << std::setw(width) << "delta sin(pi x)"
	      << std::setw(width) << "cos_pi (GCC)"
	      << std::setw(width) << "cos_pi (Boost)"
	      << std::setw(width) << "delta cos_pi"
	      << std::setw(width) << "delta cos(pi x)"
	      << std::setw(width) << "tan_pi (GCC)"
	      << std::setw(width) << "delta tan(pi x)"
	      << '\n';
    std::cout << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << std::setw(width) << "==============="
	      << '\n';
    for (int i = 0; i <= +3200; ++i)
      {
	auto x = _Tp(0.25Q * i);
	auto sin_pi_g = __gnu_cxx::sin_pi(x);
	auto sin_pi_b = beast::sin_pi(x);
	auto cos_pi_g = __gnu_cxx::cos_pi(x);
	auto cos_pi_b = beast::cos_pi(x);
	auto tan_pi_g = __gnu_cxx::tan_pi(x);
	std::cout << std::setw(width) << x
		  << std::setw(width) << sin_pi_g
		  << std::setw(width) << sin_pi_b
		  << std::setw(width) << sin_pi_g - sin_pi_b
		  << std::setw(width) << sin_pi_g - std::sin(_S_pi * x)
		  << std::setw(width) << cos_pi_g
		  << std::setw(width) << cos_pi_b
		  << std::setw(width) << cos_pi_g - cos_pi_b
		  << std::setw(width) << cos_pi_g - std::cos(_S_pi * x)
		  << std::setw(width) << tan_pi_g
		  << std::setw(width) << tan_pi_g - std::tan(_S_pi * x)
		  << '\n';
      }
    std::cout << std::endl;
  }


int
main()
{
  std::cout << "\ndouble\n=====\n\n";
  run_sin_cos_pi<double>();
}
