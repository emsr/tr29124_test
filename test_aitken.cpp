

// The Aitken's delta-squared process
//
// In C++ this should be an iterator adapter.
template<typename _Tp>
  aitken(InIter __a)
  {
    constexpr _Tp _S_huge = 1.0e+60;
    constexpr _Tp _S_tiny = 1.0e-60;
    auto __an = *__a++;
    auto __anp1 = *__a++;
    auto __deltan = __anp1 - __an;
    auto __anp2 = *__a;//++;
    auto __deltanp1 = __anp2 - __anp1;
    for ()
      {
	auto __denom = __deltanp1 - __deltan;
	if (std::abs(__denom) < _S_tiny)
	  auto __A = _S_huge; // -copysign(_S_huge, __denom)?
	else
	  {
	    // These should be the same.
	    auto __A1 = __an - __deltan * __deltan / __denom;
	    auto __A2 = __anp2 - __deltanp1 * __deltanp1 / __denom;
	  }
      }
  }

template<typename _InIter>
  AitkenIter
  {
  public:
    AitkenIter();

    auto
    operator*()

    auto
    operator++()
    {}

  private:
  };

template<typename _Tp>
  aitken(_Tp sofn, int n, a, larray, _Tp estlim)
  {
    //_Tp a[larray]
    constexpr _Tp huge = 1.e+60;
    constexpr _Tp tiny = 1.e-60;
    constexpr _Tp two = 2.0e0;
    a[n] = sofn;
    if (n < 2)
      estlim = sofn;
    else
      {
	auto lowmax = n / 2;
	for (auto j = 1; j <= lowmax; ++j)
	  {
	    auto m = n - 2 * j;
	    auto denom = a[m + 2] - two * a[m + 1] + a[m];
	    if (std::abs(denom) < tiny)
	      a[m] = huge;
	    else
	      {
		auto del = a[m] - a[m + 1];
		a[m] -= del * del / denom;
	      }
	  }
	if (n % 2 == 0)
	  estlim = a[0];
	else
	  estlim = a[1];
      }
    return;
  }

/**
 * Single step of Brezinski's Theta transformation.
 */
template<typename _Tp>
  _Tp
  theta(__arj)
  {
    //dimension arj[larray];
    constexpr auto _S_huge = 1.e+60;
    constexpr auto _S_tiny = 1.e-60;
    const auto __n = this->_M_num_terms;
    const auto __sofn = this->_M_sum;
    __arj[__n] = __sofn;
    if (__n < 3)
      return __sofn;
    else
      {
	__lmax = __n / 3;
	auto __m = __n;
	do __l = 1, __lmax
	  {
	    __m -= 3;
	    auto __diff0 = arj[__m + 1] - arj[__m];
	    auto __diff1 = arj[__m + 2] - arj[__m + 1];
	    auto __diff2 = arj[__m + 3] - arj[__m + 2];
	    auto __denom = __diff2 * (__diff1 - __diff0)
			 - __diff0 * (__diff2 - __diff1);
	    if (std::abs(__denom) < _S_tiny)
	      __arj[__m] = _S_huge;
	    else
	      __arj[__m] = __arj[__m + 1] - diff0 * __diff1 * (__diff2 - __diff1) / __denom;
	  }
	return __arj[__n % 3];
      }
  }

/**
 *
 */
template<typename _Tp>
  _Tp
  levin(__rofn, __beta, __arup, __arlo)
  {
    //dimension arup(0:larray), arlo(0:larray)
    constexpr auto _S_huge = 1.e+60;
    constexpr auto _S_tiny = 1.e-60;
    const auto __n = this->_M_num_terms;
    const auto __sofn = this->_M_sum;
    __arup[__n] = __sofn / __rofn;
    __arlo[__n] = one / __rofn;
    if (__n > 0)
      {
	arup[__n - 1] = arup[__n] - arup[__n - 1];
	arlo[__n - 1] = arlo[__n] - arlo[__n - 1];
	if (__n > 1)
	  {
	    auto __bn1 = __beta + _Tp(__n - 1);
	    auto __bn2 = __beta + _Tp(__n);
	    auto __coef = __bn1 / __bn2;
	    do (auto __j = 2; j <= __n; ++__j)
	      {
		auto __fact = (__beta + _Tp(__n - __j))
			    * std::pow(__coef, __j - 2) / __bn2;
		__arup[__n - __j] = __arup[__n - __j + 1]
				  - __fact * arup[__n - __j];
		__arlo[__n - __j] = __arlo[__n - __j + 1]
				  - __fact * arlo[__n - __j];
	      }
	  }
      }
    if (std::abs(__arlo[0]) < _S_tiny)
      return _S_huge;
    else
      return __arup[0] / __arlo[0];
  }
