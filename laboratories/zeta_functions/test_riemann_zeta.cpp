/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

#include <emsr/float128_io.h>
#include <emsr/sf_zeta.h>

#include <wrap_gsl.h>

  /**
   * @brief  Evaluate the Riemann zeta function @f$ \zeta(s) @f$
   * 	     by an alternate series for s > 0.
   *
   * The series is:
   * @f[
   * 	\zeta(s) = \frac{1}{1-2^{1-s}}
   * 		   \sum_{n=1}^{\infty} \frac{(-1)^{n-1}}{n^s}
   * @f]
   *
   * The Riemann zeta function is defined by:
   * @f[
   * 	\zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
   * @f]
   * For s < 1 use the reflection formula:
   * @f[
   * 	\zeta(s) = 2^s \pi^{s-1} \Gamma(1-s) \zeta(1-s)
   * @f]
   *
   * @todo Try vanWijnGaarden summation and maybe even Shanks avter that.
   */
  template<typename Tp>
    Tp
    riemann_zeta_alt(Tp s)
    {
      using _Val = Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto s_eps = emsr::epsilon(std::real(s));
      const unsigned int s_max_iter = 10000000;
      auto sgn = _Real{1};
      auto zeta = _Val{0};
      for (unsigned int i = 1; i < s_max_iter; ++i)
	{
	  const auto term = sgn / std::pow(_Val(i), s);
	  zeta += term;
	  sgn = -sgn;
	  if (std::abs(term) < s_eps * std::abs(zeta)
	      || std::abs(term) < s_eps
		 && std::abs(zeta) < _Real{100} * s_eps)
	    break;
	}
      zeta /= _Val{1} - std::pow(_Val{2}, _Val{1} - s);

      return zeta;
    }

  /**
   * @brief  Compute the Riemann zeta function @f$ \zeta(s) @f$
   * 	 using the product over prime factors.
   *
   * @f[
   * 	\zeta(s) = \Pi_{i=1}^\infty \frac{1}{1 - p_i^{-s}}
   * @f]
   * where @f$ {p_i} @f$ are the prime numbers.
   *
   * The Riemann zeta function is defined by:
   * @f[
   *    \renewcommand\Re{\operatorname{Re}}
   *    \renewcommand\Im{\operatorname{Im}}
   * 	\zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for \Re{s} > 1
   * @f]
   * For \Re(s) < 1 use the reflection formula:
   * @f[
   * 	\zeta(s) = (2\pi)^s \Gamma(1-s) \zeta(1-s) / \pi
   * @f]
   *
   * @param s The argument
   */
  template<typename Tp>
    Tp
    riemann_zeta_product(Tp s)
    {
      using _Val = Tp;
      using _Real = emsr::num_traits_t<_Val>;

      const auto s_eps = emsr::epsilon(std::real(s));
      constexpr unsigned long
        s_num_primes = sizeof(unsigned long) != 8 ? 256 : 256 + 48;

      auto zeta = _Val{1};
      for (unsigned long i = 0;
	   i < emsr::detail::s_num_primes; ++i)
	{
	  const auto fact = _Val{1}
			    - std::pow(_Real(emsr::prime(i)), -s);
	  zeta *= fact;
	  if (std::abs(Tp{1} - fact) < s_eps) // Assume zeta near 1.
	    break;
	}

      zeta = Tp{1} / zeta;

      return zeta;
    }

  /**
   * @brief  Return the Riemann zeta function @f$ \zeta(s) - 1 @f$
   * by summation for \Re(s) > 1.  This is a small remainder for large s.
   *
   * The Riemann zeta function is defined by:
   * @f[
   *   \zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for \Re(s) > 1
   * @f]
   *
   * @param s The argument @f$ s != 1 @f$
   */
  template<typename Tp>
    Tp
    riemann_zeta_m_1_sum(Tp s)
    {
      using _Val = Tp;
      using _Real = emsr::num_traits_t<_Val>;
      const auto s_eps = emsr::epsilon(std::real(s));
      if (emsr::fp_is_integer(s) == _Real{1})
       return emsr::quiet_NaN(std::real(s));
      else
       {
	 const auto arg = -std::log10(s_eps) / std::abs(s);
	 int k_max = arg > _Real{6} ? 1000000 : std::pow(_Real{10}, arg);
	 auto zeta_m_1 = _Val{0};
	 for (int k = k_max; k >= 2; --k)
	   {
	     auto term = std::pow(_Real(k), -s);
	     zeta_m_1 += term;
	   }
	 return zeta_m_1;
       }
    }


template<typename Tp>
  void
  plot_riemann_zeta(std::string filename, Tp proto = Tp{})
  {
    using _Val = Tp;
    using _Real = emsr::num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    const auto deg = emsr::deg_v<_Real>;

    auto data = std::ofstream(filename);

    data.precision(emsr::digits10(proto));
    data << std::showpoint << std::scientific;
    auto w = 8 + data.precision();

    using zetaT = decltype(emsr::detail::riemann_zeta(_Cmplx{}));
    std::vector<std::vector<zetaT>> sv;
    std::vector<std::vector<zetaT>> zetav;

    int i_min = -200;
    int j_min = -50;

    for (int i = i_min; i <= +50; ++i)
      {
        sv.push_back(std::vector<zetaT>{});
	zetav.push_back(std::vector<zetaT>{});
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = _Cmplx(0.10L * i, 0.10L * j);
	    sv.back().push_back(s);
	    zetav.back().push_back(emsr::detail::riemann_zeta(s));
	  }
      }

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = sv[i - i_min][j - j_min];
	    auto zeta = zetav[i - i_min][j - j_min];
	    data << ' ' << std::setw(w) << std::real(s)
		 << ' ' << std::setw(w) << std::imag(s)
		 << ' ' << std::setw(w) << std::real(zeta)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = sv[i - i_min][j - j_min];
	    auto zeta = zetav[i - i_min][j - j_min];
	    data << ' ' << std::setw(w) << std::real(s)
		 << ' ' << std::setw(w) << std::imag(s)
		 << ' ' << std::setw(w) << std::imag(zeta)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = sv[i - i_min][j - j_min];
	    auto zeta = zetav[i - i_min][j - j_min];
	    data << ' ' << std::setw(w) << std::real(s)
		 << ' ' << std::setw(w) << std::imag(s)
		 << ' ' << std::setw(w) << std::abs(zeta)
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';

    for (int i = i_min; i <= +50; ++i)
      {
	for (int j = j_min; j <= +50; ++j)
	  {
	    auto s = sv[i - i_min][j - j_min];
	    auto zeta = zetav[i - i_min][j - j_min];
	    data << ' ' << std::setw(w) << std::real(s)
		 << ' ' << std::setw(w) << std::imag(s)
		 << ' ' << std::setw(w) << deg * std::arg(zeta) 
		 << '\n';
	  }
	data << '\n';
      }
    data << '\n';
  }

template<typename Tp>
  void
  test_riemann_zeta_real(Tp proto = Tp{})
  {
    using _Val = Tp;
    using _Real = emsr::num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    int i_min = -250;

    std::cout << '\n'
	      << ' ' << std::setw(w) << "s"
	      << ' ' << std::setw(4 + 2 * w) << "zetac = zeta (cmplx)"
	      << ' ' << std::setw(w) << "zeta (real)"
	      << ' ' << std::setw(w) << "|zetac - zeta|"
	      << '\n';
    for (int i = i_min; i <= +250; ++i)
      {
        auto sc = _Cmplx(0.10L * i, 0.0L);
	auto zetac = emsr::detail::riemann_zeta(sc);
        auto s = _Real(0.10L * i);
	auto zeta = emsr::detail::riemann_zeta(s);
	std::cout << ' ' << std::setw(w) << s
		  << ' ' << std::setw(4 + 2 * w) << zetac
		  << ' ' << std::setw(w) << zeta
		  << ' ' << std::setw(w) << std::abs(zetac - zeta)
		  << '\n';
      }
  }

template<typename Tp>
  void
  test_riemann_zeta(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    int i_min = -500;

    std::cout << '\n'
	      << ' ' << std::setw(w) << "s"
	      << ' ' << std::setw(w) << "zeta"
	      << ' ' << std::setw(w) << "zeta_GSL"
	      << ' ' << std::setw(w) << "|zeta - zeta_GSL|"
	      << '\n';
    for (int i = i_min; i <= +500; ++i)
      {
        auto s = Tp(0.05L * i);
	if (s == Tp{1})
	  {
	    std::cout << ' ' << std::setw(w) << s
		      << ' ' << std::setw(w) << "nan"
		      << ' ' << std::setw(w) << "nan"
		      << ' ' << std::setw(w) << "nan"
		      << '\n';
	    continue;
	  }
	auto zeta = emsr::riemann_zeta(s);
	auto zeta_gsl = gsl::riemann_zeta(s);
	std::cout << ' ' << std::setw(w) << s
		  << ' ' << std::setw(w) << zeta
		  << ' ' << std::setw(w) << zeta_gsl
		  << ' ' << std::setw(w) << zeta - zeta_gsl
		  << '\n';
      }
  }

template<typename Tp>
  void
  test_nontrivial_zeros(Tp proto = Tp{})
  {
    using namespace std::complex_literals;

    using _Val = Tp;
    using _Real = emsr::num_traits_t<_Val>;
    using _Cmplx = std::complex<_Real>;

    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    int i_min = -500;

    std::cout << '\n'
	      << ' ' << std::setw(w) << "s"
	      << ' ' << std::setw(w) << "zeta"
	      << '\n';
    for (int i = i_min; i <= +500; ++i)
      {
        auto s = _Cmplx(0.5L, 0.05L * i);
	auto zeta = emsr::detail::riemann_zeta(s);
	std::cout << ' ' << std::setw(w) << std::imag(s)
		  << ' ' << std::setw(w) << std::real(zeta)
		  << ' ' << std::setw(w) << std::imag(zeta)
		  << ' ' << std::setw(w) << std::abs(zeta)
		  << '\n';
      }

    _Cmplx
    zeros[10]
    {
      0.5l + 14.134725141734693790457251983562470270784257115699il,
      0.5l + 21.022039638771554992628479593896902777334340524903il,
      0.5l + 25.010857580145688763213790992562821818659549672558il,
      0.5l + 30.424876125859513210311897530584091320181560023715il,
      0.5l + 32.935061587739189690662368964074903488812715603517il,
      0.5l + 37.586178158825671257217763480705332821405597350831il,
      0.5l + 40.918719012147495187398126914633254395726165962777il,
      0.5l + 43.327073280914999519496122165406805782645668371837il,
      0.5l + 48.005150881167159727942472749427516041686844001144il,
      0.5l + 49.773832477672302181916784678563724057723178299677il
    };

    std::cout << '\n'
	      << ' ' << std::setw(4 + 2 * w) << "s"
	      << ' ' << std::setw(4 + 2 * w) << "zeta(s)"
	      << ' ' << std::setw(w) << "|zeta(s)|"
	      << '\n';
    for (auto s : zeros)
      {
	auto zeta = emsr::detail::riemann_zeta(s);
	std::cout << ' ' << std::setw(4 + 2 * w) << s
		  << ' ' << std::setw(4 + 2 * w) << zeta
		  << ' ' << std::setw(w) << std::abs(zeta)
		  << '\n';
      }
  }

int
main(int n_app_args, char** arg)
{
  using namespace std::complex_literals;

  // These barf on Cygwin because the literals get turned to complex__ long double!
  //auto zetam = emsr::detail::riemann_zeta(0.01l - 1.0il);
  //std::cout << "zeta(" << 0.01l - 1.0il << ") = " << zetam << '\n';
  //auto zetap = emsr::detail::riemann_zeta(0.01l + 1.0il);
  //std::cout << "zeta(" << 0.01l + 1.0il << ") = " << zetap << '\n';

  test_riemann_zeta(1.0);

  test_riemann_zeta_real<long double>();

  test_nontrivial_zeros<long double>();

  std::cout << "\n\nRiemann zeta\n\n" << std::flush;

  std::string plot_data_dir = ".";
  if (n_app_args > 1)
    plot_data_dir = arg[1];

  std::cout << "\nriemann_zeta<float>\n";
  plot_riemann_zeta<float>(plot_data_dir + '/' + "riemann_zeta_float.txt");

  std::cout << "\nriemann_zeta<double>\n";
  plot_riemann_zeta<double>(plot_data_dir + '/' + "riemann_zeta_double.txt");

  std::cout << "\nriemann_zeta<long double>\n";
  plot_riemann_zeta<long double>(plot_data_dir + '/' + "riemann_zeta_long_double.txt");

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //std::cout << "\nriemann_zeta<__float128>\n";
  //plot_riemann_zeta<__float128>(plot_data_dir + '/' + "riemann_zetafloat128.txt");
#endif
}
