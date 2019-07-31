/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>

template<typename _Tp>
  void
  test()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    auto w = 8 + std::cout.precision();

    const auto _S_ln2 = __gnu_cxx::numbers::__ln_2_v<_Tp>;

    //std::cout << std::setfill('0');

    std::cout << "ln(2) = " << std::log(_Tp{2}) << '\n';
    std::cout << "ln(2) = " << _S_ln2 << '\n';

    for (auto __l = 1u; __l < std::__detail::_S_num_factorials<_Tp>; __l += 2)
      {
	std::cout << '\n';
	  {
	    const auto __lm = __l - 1;
	    const auto __lmfact = std::__detail::__factorial<_Tp>(__lm);
	    const auto __mm = __lm / 2;
	    const auto __mmfact = std::__detail::__factorial<_Tp>(__mm);
	    auto __Plm1 = ((__mm & 1) ? -1 : 1) * __lmfact / __mmfact / __mmfact
			  / std::pow(_Tp{2}, __lm);
	    auto __Ppl = __l * __Plm1;
	    auto __weight = _Tp{2} / __Ppl / __Ppl;
	    std::cout << ' ' << std::setw(w) << __Plm1 << ' ' << std::setw(w) << __Ppl << ' ' << std::setw(w) << __weight << '\n';
	  }

	  {
	    const auto __lm = __l - 1;
	    const auto __lmfact = std::__detail::__log_factorial<_Tp>(__lm);
	    const auto __mm = __lm / 2;
	    const auto __mmfact = std::__detail::__log_factorial<_Tp>(__mm);
	    auto __Plm1 = ((__mm & 1) ? -1 : 1)
			* std::exp(__lmfact - 2 * __mmfact - __lm * _S_ln2);

	    auto __Ppl = __l * __Plm1;
	    auto __weight = _Tp{2} / __Ppl / __Ppl;
	    std::cout << ' ' << std::setw(w) << __Plm1 << ' ' << std::setw(w) << __Ppl << ' ' << std::setw(w) << __weight << '\n';
	  }

	  {
	    const auto __lm = __l - 1;
	    const auto __mm = __lm / 2;
	    auto _Am = _Tp{1};
	    for (auto __m = 1u; __m <= __mm; ++__m)
	      _Am *= -_Tp(2 * __m - 1) / _Tp(2 * __m);
	    auto __Plm1 = _Am;
	    auto __Ppl = __l * __Plm1;
	    auto __weight = _Tp{2} / __Ppl / __Ppl;
	    std::cout << ' ' << std::setw(w) << __Plm1 << ' ' << std::setw(w) << __Ppl << ' ' << std::setw(w) << __weight << '\n';
	  }
      }
  }

int
main()
{
  test<double>();
}
