

  /**
   *  @brief Evaluate polynomial based on the confluent hypergeometric
   *         representation.
   *  @f[
   *    L_n^\alpha(x) = (\alpha + 1)_n / n! _1F_1(-n, \alpha + 1, x)
   *  @f]
   *  Assumes n > 0 and @f$ \alpha @f$ != negative integer greater than -n.
   */
  template<typename _Tpa, class _Tp>
  _Tp
    __laguerre_n_cp(unsigned int __n, _Tpa __alpha, _Tp __x)
    {
      constexpr auto __max = _Tp{0.9L} * std::numeric_limits<_Tp>::max();
      const auto __lnfact = __log_gamma(__n + _Tp{1});
      const auto __lg1 = __log_gamma(__alpha + _Tp{1} + __n);
      const auto __s1 = __log_gamma_sign(__alpha + _Tp{1} + __n);
      const auto __lg2 = __log_gamma(__alpha + _Tp{1});
      const auto __s2 = __log_gamma_sign(__alpha + _Tp{1});
      const auto lnpre = (__lg1 - __lg2) - __lnfact;

      auto __poly_1F1 = _Tp{1};
      for (auto __k = int(__n) - 1; __k >= 0; --__k)
	{
	  auto __t = (-int(__n) + __k) * (__x / _Tp(__k + 1))
		   / (__alpha + _Tp{1} + __k);
	  auto __r = __t + _Tp{1} / __poly_1F1;
	  if (__r > __max / __poly_1F1) // Internal error only, catch in the main routine.
            std::__throw_runtime_error(__N("__laguerre_n_cp: series failed));
	  else // Collect the Horner terms.
            __poly_1F1  = _Tp{1} + __t * __poly_1F1;
	}

      return std::exp(__lnpre) * __poly_1F1;
    }
