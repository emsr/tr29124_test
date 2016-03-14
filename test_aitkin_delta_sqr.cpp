

// The Aitken's delta-squared process
//
// In C++ this should be an iterator adapter.
template<typename _Tp>
  aitken(InIter __a)
  {
    auto __an = *__a++;
    auto __anp1 = *__a++;
    auto __deltan = __anp1 - __an;
    auto __anp2 = *__a++;
    auto __deltanp1 = __anp2 - __anp1;
    for ()
      {
	auto __A1 = __an - __deltan * __deltan / (__deltanp1 - __deltan);
	auto __A2 = __anp2 - __deltanp1 * __deltanp1 / (__deltanp1 - __deltan);
      }
  }
