

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
