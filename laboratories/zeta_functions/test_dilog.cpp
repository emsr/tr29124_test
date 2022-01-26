/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>

#ifdef NO_DILOG

namespace emsr
{
namespace detail
{

  template<typename _Tp>
    _Tp
    dilog(_Tp x)
    {
      static constexpr unsigned long long s_maxit = 100000ULL;
      static const auto s_eps = 10 * emsr::epsilon(x);
      static const auto s_pipio6 = emsr::pi_sqr_div_6_v<_Tp>;
      if (std::isnan(x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (x > _Tp(+1))
	throw std::range_error(__N("dilog: argument greater than one"));
      else if (x < _Tp(-1))
	{
	  auto lnfact = std::log(_Tp(1) - x);
	  return -dilog(_Tp(1) - _Tp(1) / (_Tp(1) - x))
		 - _Tp(0.5L) * lnfact * lnfact;
	}
      else if (x == _Tp(1))
	return s_pipio6;
      else if (x == -_Tp(1))
	return -_Tp(0.5L) * s_pipio6;
      else if (x > _Tp(0.5))
	return s_pipio6 - std::log(x) * std::log(_Tp(1) - x)
	     - dilog(_Tp(1) - x);
      else if (x < -_Tp(0.5))
	return -_Tp(0.5L) * s_pipio6 - std::log(_Tp(1) + x) * std::log(-x)
	     + dilog(_Tp(1) + x) - dilog(_Tp(1) - x * x);
      else
	{
	  _Tp sum = 0;
	  _Tp fact = 1;
	  for (auto i = 1ULL; i < s_maxit; ++i)
	    {
	      fact *= x;
	      auto term = fact / (i * i);
	      sum += term;
	      if (std::abs(term) < s_eps)
		{
		  std::cout << i << " - ";
		  break;
		}
	      if (i + 1 == s_maxit)
		throw std::runtime_error("dilog: sum failed");
	    }
	  return sum;
	}
    }

} // namespace detail
} // namespace emsr

namespace emsr
{

  //  Dilogarithm functions

  inline float
  dilogf(float x)
  { return emsr::detail::dilog<float>(x); }

  inline long double
  dilogl(long double x)
  { return emsr::detail::dilog<long double>(x); }

  template<typename _Tp>
    inline typename emsr::promote<_Tp>::type
    dilog(_Tp x)
    {
      typedef typename emsr::promote<_Tp>::type type;
      return emsr::detail::dilog<type>(x);
    }

}

#endif // NO_DILOG

int
main()
{
  std::cout.precision(15);
  std::cout << '\n';
  for (auto i = -4000; i <= +1000; ++i)
    std::cout << "dilog(" << i * 0.001 << ") = " << emsr::dilog(i * 0.001) << '\n';
}
