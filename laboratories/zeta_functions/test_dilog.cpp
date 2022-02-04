/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>

#include <emsr/sf_zeta.h>

#ifdef NO_DILOG

namespace emsr
{
namespace detail
{

  template<typename Tp>
    Tp
    dilog(Tp x)
    {
      static constexpr unsigned long long s_maxit = 100000ULL;
      static const auto s_eps = 10 * emsr::epsilon(x);
      static const auto s_pipio6 = emsr::pi_sqr_div_6_v<Tp>;
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else if (x > Tp(+1))
	throw std::range_error(__N("dilog: argument greater than one"));
      else if (x < Tp(-1))
	{
	  auto lnfact = std::log(Tp(1) - x);
	  return -dilog(Tp(1) - Tp(1) / (Tp(1) - x))
		 - Tp(0.5L) * lnfact * lnfact;
	}
      else if (x == Tp(1))
	return s_pipio6;
      else if (x == -Tp(1))
	return -Tp(0.5L) * s_pipio6;
      else if (x > Tp(0.5))
	return s_pipio6 - std::log(x) * std::log(Tp(1) - x)
	     - dilog(Tp(1) - x);
      else if (x < -Tp(0.5))
	return -Tp(0.5L) * s_pipio6 - std::log(Tp(1) + x) * std::log(-x)
	     + dilog(Tp(1) + x) - dilog(Tp(1) - x * x);
      else
	{
	  Tp sum = 0;
	  Tp fact = 1;
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

  template<typename Tp>
    inline typename emsr::promote<Tp>::type
    dilog(Tp x)
    {
      typedef typename emsr::promote<Tp>::type type;
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
