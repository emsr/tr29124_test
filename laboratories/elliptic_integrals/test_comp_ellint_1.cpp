/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <utility>

#include <emsr/float128_io.h>

//  Use AGM to do an ab initio calculation of K(k).
template<typename _Tp>
  _Tp
  comp_ellint_1_agm(_Tp k)
  {
    constexpr auto s_max_iter = 100;
    auto am = (_Tp{1} + k) / _Tp{2};
    auto gm = std::sqrt(k);
    for (int i = 0; i < s_max_iter; ++i)
      {
	gm = std::sqrt(gm * std::exchange(am, (am + gm) / _Tp{2}));
	if (std::abs(am - gm) < std::numeric_limits<_Tp>::epsilon())
	  break;
      }
    return gm;
  }

/**
 *  Use MacLaurin series to calculate the elliptic nome
 *  given the , k.
 */
template<typename _Tp>
  _Tp
  ellint_nome_series(_Tp k)
  {
    auto m = k * k; 
    return m * ((_Tp{1} / _Tp{16})
	 + m * ((_Tp{1} / _Tp{32})
	 + m * ((_Tp{21} / _Tp{1024})
	 + m * ((_Tp{31} / _Tp{2048})
	 + m * (_Tp{6257} / _Tp{524288})))));
  }

/**
 *  Use the arithmetic-geometric mean to calculate the elliptic nome
 *  given the , k.
 */
template<typename _Tp>
  _Tp
  ellint_nome_agm(_Tp k)
  {
    const auto s_pi = emsr::pi_v<_Tp>;
    auto kp = std::sqrt((_Tp{1} - k) * (_Tp{1} + k));
    auto K = comp_ellint_1_agm(k);
    auto Kp = comp_ellint_1_agm(kp);
    return std::exp(-s_pi * Kp / K);
  }

/**
 *  Return the elliptic nome given the , k.
 */
template<typename _Tp>
  _Tp
  ellint_nome(_Tp k)
  {
    const auto s_eps = std::numeric_limits<_Tp>::epsilon();
    if (k < std::pow(_Tp{67} * s_eps, _Tp{1}/_Tp{8}))
      return ellint_nome_series(k);
    else
      return ellint_nome_agm(k);
  }


template<typename _Tp>
  void
  test_K(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));

    auto width = 6 + std::cout.precision();
    const auto del = _Tp{1} / _Tp{1000};
    for (int i = 0; i <= 1000; ++i)
      {
	auto k = i * del;
	std::cout << ' ' << std::setw(width) << k
		  << ' ' << std::setw(width) << comp_ellint_1_agm(k)
		  << ' ' << std::setw(width) << ellint_nome(k) << '\n';
      }
  }


int
main()
{
  std::cout << "\n\nfloat\n";
  test_K<float>();

  std::cout << "\n\ndouble\n";
  test_K<double>();

  std::cout << "\n\nlong double\n";
  test_K<long double>();

  //std::cout << "\n\n__float128\n";
  //test_K<__float128>();
}
