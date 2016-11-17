#ifndef STEEDCONTINUEDFRACTION_TCC
#define STEEDCONTINUEDFRACTION_TCC 1

/**
 * A Steed continued fraction evaluator.
 * This class required three functors:
 *   A partial numerator function @f$ a_k(x) @f$
 *   A partial denominator function @f$ b_k(x) @f$
 *   A tail function @f$ w_n(x) @f$
 */
template<typename _Tp, typename _AFun, typename _BFun, typename _TailFun>
  class _SteedContinuedFraction
  {
  private:

    _AFun _M_num;
    _BFun _M_den;
    _TailFun _M_tail;

  public:

    using _ARet = decltype(_M_num(0ull, _Tp{}));
    using _BRet = decltype(_M_den(0ull, _Tp{}));
    using _Ret = decltype(_M_num(0ull, _Tp{}) / _M_den(0ull, _Tp{}));

    constexpr _SteedContinuedFraction(_AFun __a, _BFun __b, _TailFun __w)
    : _M_num(__a),
      _M_den(__b),
      _M_tail(__w)
    { }

    _Ret
    operator()(_Tp __x) const
    {
      const auto _S_fp_min = 1000 * __gnu_cxx::__min(__x);
      const auto _S_eps = _Tp{0.125L} * __gnu_cxx::__epsilon(__x);
      constexpr std::size_t _S_max_iter = 1000;

      auto __b = _M_den(0, __x);
      _Ret __C(__b);
      __b = _M_den(1, __x);
      _Ret __D = _Tp{1} / __b;
      auto __a = _M_num(1, __x);
      _Ret __delta_C = __a * __D;
      __C += __delta_C;
      std::size_t __i = 2;
      while (true)
	{
	  __a = _M_num(__i, __x);
	  __b = _M_den(__i, __x);
	  __D = __a * __D + __b;
	  if (std::abs(__D) < _S_fp_min)
	    __D = _S_fp_min;
	  __D = _Ret{1} / __D;
	  __delta_C *= __b * __D - _Ret{1};
	  __C += __delta_C;
	  if (std::abs(__delta_C) < _S_eps * std::abs(__C))
	    return __C;
	  if (__i > _S_max_iter)
	    throw std::runtime_error("_SteedContinuedFraction: " "continued fraction evaluation failed");
	    //std::__throw_runtime_error(__N"_SteedContinuedFraction: "
		//		   "continued fraction evaluation failed"));
	  ++__i;
	}
    }
  };

#endif // STEEDCONTINUEDFRACTION_TCC
