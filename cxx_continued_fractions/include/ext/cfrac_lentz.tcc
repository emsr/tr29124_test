#ifndef CFRAC_LENTZ_TCC
#define CFRAC_LENTZ_TCC 1

#include <stdexcept>
#include <bits/numeric_limits.h>

/**
 * A modified Lentz continued fraction evaluator.
 * This class requires three functors:
 *   A partial numerator function @f$ a_k(x) @f$
 *   A partial denominator function @f$ b_k(x) @f$
 *   A tail function @f$ w_n(x) @f$
 */
template<typename _Tp, typename _AFun, typename _BFun, typename _TailFun>
  class _LentzContinuedFraction
  {
  private:

    _AFun _M_num;
    _BFun _M_den;
    _TailFun _M_tail;

  public:

    using _ARet = decltype(_M_num(0ull, _Tp{}));
    using _BRet = decltype(_M_den(0ull, _Tp{}));
    using _Ret = decltype(_M_num(0ull, _Tp{}) / _M_den(0ull, _Tp{}));
    using _Val = __gnu_cxx::fp_promote_t<_Tp, _ARet, _BRet>;
    using _Real = std::__detail::__num_traits_t<_Val>;

    const _Real _S_fp_min = _Real{1000} * std::numeric_limits<_Real>::min();

    _Real _M_rel_error = std::numeric_limits<_Real>::epsilon();
    std::size_t _M_max_iter = 1000;

    constexpr _LentzContinuedFraction(_AFun __a, _BFun __b, _TailFun __w)
    : _M_num(__a),
      _M_den(__b),
      _M_tail(__w)
    { }

    _Ret
    operator()(_Tp __x) const
    {
      auto __b = _M_den(0, __x);
      _Ret __C(__b);
      if (std::abs(__C) < _S_fp_min)
	__C = _S_fp_min;
      auto __D = _Ret{0};
      auto __E = __C;
      std::size_t __k = 1;
      while (true)
	{
	  auto __a = _M_num(__k, __x);
	  __b = _M_den(__k, __x);
	  __D = __a * __D + __b;
	  if (std::abs(__D) < _S_fp_min)
	    __D = _S_fp_min;
	  __D = _Ret{1} / __D;
	  __E = __b + __a / __E;
	  if (std::abs(__E) < _S_fp_min)
	    __E = _S_fp_min;
	  auto __H = __E * __D;
	  __C *= __H;
	  if (std::abs(__H - _Ret{1}) < _M_rel_error)
	    return __C;
	  ++__k;
	  if (__k > _M_max_iter)
	    throw std::runtime_error("_LentzContinuedFraction: "
				     "continued fraction evaluation failed");
	    //std::__throw_runtime_error(__N"_LentzContinuedFraction: "
		//		   "continued fraction evaluation failed"));
	}
    }
  };

#endif // CFRAC_LENTZ_TCC
