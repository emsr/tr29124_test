#ifndef CFRAC_FORWARD_TCC
#define CFRAC_FORWARD_TCC 1

#include <stdexcept>

#include <emsr/fp_type_util.h>
#include <emsr/complex_util.h>
#include <emsr/numeric_limits.h>

namespace emsr
{

/**
 * You need a functors returning the numerator and denominator a, b
 * an iterator pair for a and an iterator for b
 *
 */
template<typename Tp, typename AFun, typename BFun, typename TailFun>
  class ForwardContinuedFraction
  {
  private:

    AFun m_num;
    BFun m_den;
    TailFun m_tail;

  public:

    using ARet = decltype(m_num(0ull, Tp{}));
    using BRet = decltype(m_den(0ull, Tp{}));
    using Ret = decltype(ARet(0, Tp{}) / BRet(0, Tp{}));
    using Val = emsr::fp_promote_t<Tp, ARet, BRet>;
    using Real = emsr::num_traits_t<Val>;

    Real m_rel_error = std::numeric_limits<Real>::epsilon();
    int m_max_iter = 1000;

    ForwardContinuedFraction(AFun a, BFun b, TailFun w)
    : m_num(a),
      m_den(b),
      m_tail(w)
    { }

    Ret
    operator()(Tp x) const
    {
      auto b = m_den(0, x);
      auto Den2 = b;
      auto a = m_num(1, x);
      auto Num2 = a;
      b = m_den(2, x);
      auto Num1 = b * Num2;
      a = m_num(2, x);
      auto Den1 = a + b * Den2;
      std::size_t k = 2;
      while (true)
	{
	  a = m_num(k, x);
	  b = m_den(k, x);
	  auto Num = b * Num1 + a * Num2;
	  auto Den = b * Den1 + a * Den2;
	  Num2 = Num1;
	  Den2 = Den1;
	  Num1 = Num;
	  Den1 = Den;
	  ++k;
	  if (k > m_max_iter)
	    throw std::runtime_error("ForwardContinuedFraction: "
				     "continued fraction evaluation failed");
	}
    }
  };

} // namespace emsr

#endif // CFRAC_FORWARD_TCC
