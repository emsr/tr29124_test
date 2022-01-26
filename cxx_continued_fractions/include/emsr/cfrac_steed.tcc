#ifndef CFRAC_STEED_TCC
#define CFRAC_STEED_TCC 1

#include <emsr/fp_type_util.h>
#include <emsr/complex_util.h>

namespace emsr
{

/**
 * A Steed continued fraction evaluator.
 * This class requires three functors:
 *   A partial numerator function @f$ a_k(x) @f$
 *   A partial denominator function @f$ b_k(x) @f$
 *   A tail function @f$ w_n(x) @f$
 */
template<typename Tp, typename AFun, typename BFun, typename TailFun>
  class SteedContinuedFraction
  {
  private:

    AFun m_num;
    BFun m_den;
    TailFun m_tail;

  public:

    using ARet = decltype(m_num(0ull, Tp{}));
    using BRet = decltype(m_den(0ull, Tp{}));
    using Ret = decltype(m_num(0ull, Tp{}) / m_den(0ull, Tp{}));
    using Val = emsr::fp_promote_t<Tp, ARet, BRet>;
    using Real = emsr::num_traits_t<Val>;

    const Real s_fp_min = Real{1000} * std::numeric_limits<Real>::min();

    Real m_rel_error = std::numeric_limits<Real>::epsilon();
    std::size_t m_max_iter = 1000;

    constexpr SteedContinuedFraction(AFun a, BFun b, TailFun w)
    : m_num(a),
      m_den(b),
      m_tail(w)
    { }

    Ret
    operator()(Tp x) const
    {
      auto b = m_den(0, x);
      Ret C(b);
      b = m_den(1, x);
      Ret D = Tp{1} / b;
      auto a = m_num(1, x);
      Ret delta_C = a * D;
      C += delta_C;
      std::size_t k = 2;
      while (true)
	{
	  a = m_num(k, x);
	  b = m_den(k, x);
	  D = a * D + b;
	  if (std::abs(D) < s_fp_min)
	    D = s_fp_min;
	  D = Ret{1} / D;
	  delta_C *= b * D - Ret{1};
	  C += delta_C;
	  if (std::abs(delta_C) < m_rel_error * std::abs(C))
	    return C;
	  ++k;
	  if (k > m_max_iter)
	    throw std::runtime_error("SteedContinuedFraction: continued fraction evaluation failed");
	    //std::throw_runtime_error("SteedContinuedFraction: continued fraction evaluation failed");
	}
    }
  };

} // namespace emsr

#endif // CFRAC_STEED_TCC
