/**
 *
 */

#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/math_constants.h>
#include <emsr/sf_gamma.h> / factorial

template<typename _Tp>
  void
  test()
  {
    std::cout.precision(emsr::digits10<_Tp>());
    auto w = 8 + std::cout.precision();

    const auto _S_ln2 = emsr::ln2_v<_Tp>;

    //std::cout << std::setfill('0');

    std::cout << "ln(2) = " << std::log(_Tp{2}) << '\n';
    std::cout << "ln(2) = " << _S_ln2 << '\n';

    for (auto l = 1u; l < emsr::detail::s_num_factorials<_Tp>; l += 2)
      {
	std::cout << '\n';
	  {
	    const auto lm = l - 1;
	    const auto lmfact = emsr::detail::factorial<_Tp>(lm);
	    const auto mm = lm / 2;
	    const auto mmfact = emsr::detail::factorial<_Tp>(mm);
	    auto Plm1 = ((mm & 1) ? -1 : 1) * lmfact / mmfact / mmfact
			  / std::pow(_Tp{2}, lm);
	    auto Ppl = l * Plm1;
	    auto weight = _Tp{2} / Ppl / Ppl;
	    std::cout << ' ' << std::setw(w) << Plm1 << ' ' << std::setw(w) << Ppl << ' ' << std::setw(w) << weight << '\n';
	  }

	  {
	    const auto lm = l - 1;
	    const auto lmfact = emsr::detail::log_factorial<_Tp>(lm);
	    const auto mm = lm / 2;
	    const auto mmfact = emsr::detail::log_factorial<_Tp>(mm);
	    auto Plm1 = ((mm & 1) ? -1 : 1)
			* std::exp(lmfact - 2 * mmfact - lm * _S_ln2);

	    auto Ppl = l * Plm1;
	    auto weight = _Tp{2} / Ppl / Ppl;
	    std::cout << ' ' << std::setw(w) << Plm1 << ' ' << std::setw(w) << Ppl << ' ' << std::setw(w) << weight << '\n';
	  }

	  {
	    const auto lm = l - 1;
	    const auto mm = lm / 2;
	    auto _Am = _Tp{1};
	    for (auto m = 1u; m <= mm; ++m)
	      _Am *= -_Tp(2 * m - 1) / _Tp(2 * m);
	    auto Plm1 = _Am;
	    auto Ppl = l * Plm1;
	    auto weight = _Tp{2} / Ppl / Ppl;
	    std::cout << ' ' << std::setw(w) << Plm1 << ' ' << std::setw(w) << Ppl << ' ' << std::setw(w) << weight << '\n';
	  }
      }
  }

int
main()
{
  test<double>();
}
