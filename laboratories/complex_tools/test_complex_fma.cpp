/**
 *
 */

#include <cmath>
#include <complex>
#include <limits>

  template<typename _Tp>
    std::complex<_Tp>
    fma_naive(const std::complex<_Tp>& a, const std::complex<_Tp>& b,
	      const std::complex<_Tp>& c)
    { return a * b + c; }

  /*
   * This is bogus because it has normal adds in it.
   */
  template<typename _Tp>
    std::complex<_Tp>
    fma_sym(const std::complex<_Tp>& a, const std::complex<_Tp>& b,
	    const std::complex<_Tp>& c)
    {
      const auto [ar, ai] = reinterpret_cast<const _Tp(&)[2]>(a);
      const auto [br, bi] = reinterpret_cast<const _Tp(&)[2]>(b);
      const auto [cr, ci] = reinterpret_cast<const _Tp(&)[2]>(c);
      // Try to symmetrize the treatment of Im(c).  We need two terms anyway.
      return std::complex<_Tp>(
	std::fma(-ai, bi, std::fma(ar, br, cr)),
	std::fma(ar, bi, _Tp{0.5L} * ci)
	 + std::fma(ai, br, ci - _Tp{0.5L} * ci));
    }

  template<typename _Tp>
    std::complex<_Tp>
    fma(const std::complex<_Tp>& a, const std::complex<_Tp>& b,
	const std::complex<_Tp>& c)
    {
      const auto [ar, ai] = reinterpret_cast<const _Tp(&)[2]>(a);
      const auto [br, bi] = reinterpret_cast<const _Tp(&)[2]>(b);
      const auto [cr, ci] = reinterpret_cast<const _Tp(&)[2]>(c);
      return std::complex<_Tp>(
	std::fma(ar, br, -std::fma(ai, bi, -cr)),
	std::fma(ar, bi, std::fma(ai, br, ci)));
    }

template<typename _Tp>
  void
  test_fma()
  {
    using cmplx = std::complex<_Tp>;
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    const auto wr = std::cout.precision() + 8;
    const auto wc = 2 * wr + 4;

    std::vector<cmplx> aa;
    for (int i = -3; i <= +3; ++i)
      for (int j = -3; j <= +3; ++j)
        aa.push_back(cmplx(i / _Tp{3}, j / _Tp{3}));

    std::vector<cmplx> bb;
    for (int i = -5; i <= +5; ++i)
      for (int j = -5; j <= +5; ++j)
        bb.push_back(cmplx(i / _Tp{5}, j / _Tp{5}));

    for (auto&& a : aa)
      for (auto&& b : bb)
	for (int i = -10; i <= +10; ++i)
	  {
	    if (i == -10)
	      std::cout << '\n';
	    for (int j = -10; j <= +10; ++j)
	      {
		const auto c = cmplx(_Tp{0.5L} * i, _Tp{0.5L} * j);
		const auto naive = fma_naive(a, b, c);
		//const auto symm = fma_sym(a, b, c);
		const auto straight = fma(a, b, c);
		std::cout << ' ' << std::setw(wc) << a
			  << ' ' << std::setw(wc) << b
			  << ' ' << std::setw(wc) << c
			  << ' ' << std::setw(wc) << naive
		//	  << ' ' << std::setw(wc) << symm
			  << ' ' << std::setw(wc) << straight
			  << ' ' << std::setw(wr) << std::abs(naive - straight)
		//	  << ' ' << std::setw(wr) << std::abs(symm - straight)
			  << '\n';
	      }
	  }
  }

int
main()
{
  //test_fma<float>();

  test_fma<double>();

  test_fma<long double>();

  //test_fma<>();

  return 0;
}
