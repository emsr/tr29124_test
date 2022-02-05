/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <utility>

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>
#include <emsr/float128_io.h>

//  Use AGM to do an ab initio calculation of K(k).
template<typename Tp>
  Tp
  comp_ellint_1_agm(Tp k)
  {
    constexpr auto s_max_iter = 100;
    auto am = (Tp{1} + k) / Tp{2};
    auto gm = std::sqrt(k);
    for (int i = 0; i < s_max_iter; ++i)
      {
	gm = std::sqrt(gm * std::exchange(am, (am + gm) / Tp{2}));
	if (std::abs(am - gm) < std::numeric_limits<Tp>::epsilon())
	  break;
      }
    return gm;
  }

/**
 *  Use MacLaurin series to calculate the elliptic nome
 *  given the , k.
 */
template<typename Tp>
  Tp
  ellint_nome_series(Tp k)
  {
    auto m = k * k; 
    return m * ((Tp{1} / Tp{16})
	 + m * ((Tp{1} / Tp{32})
	 + m * ((Tp{21} / Tp{1024})
	 + m * ((Tp{31} / Tp{2048})
	 + m * (Tp{6257} / Tp{524288})))));
  }

/**
 *  Use the arithmetic-geometric mean to calculate the elliptic nome
 *  given the , k.
 */
template<typename Tp>
  Tp
  ellint_nome_agm(Tp k)
  {
    const auto s_pi = emsr::pi_v<Tp>;
    auto kp = std::sqrt((Tp{1} - k) * (Tp{1} + k));
    auto K = comp_ellint_1_agm(k);
    auto Kp = comp_ellint_1_agm(kp);
    return std::exp(-s_pi * Kp / K);
  }

/**
 *  Return the elliptic nome given the , k.
 */
template<typename Tp>
  Tp
  ellint_nome(Tp k)
  {
    const auto s_eps = std::numeric_limits<Tp>::epsilon();
    if (k < std::pow(Tp{67} * s_eps, Tp{1}/Tp{8}))
      return ellint_nome_series(k);
    else
      return ellint_nome_agm(k);
  }


template<typename Tp>
  void
  test_K(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));

    auto width = 6 + std::cout.precision();
    const auto del = Tp{1} / Tp{1000};
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
