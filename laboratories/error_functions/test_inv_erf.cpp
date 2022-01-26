/**
 *
 */

// AAOF pp. 408-409.

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>

#include <emsr/float128_io.h>
namespace detail
{

  template<typename _Tp>
    _Tp
    erfc_scaled(_Tp x)
    { return std::exp(x * x) * std::erfc(x); }

  /**
   * Return the inverse error function by recursion:
   * @f[
   *
   * @f]
   *
   * @param[in] p The argument between -1 and 1
   * @param[in] x The initial x-value guess.
   */
  template<typename _Tp>
    _Tp
    erf_inv_recur(_Tp p, _Tp x = _Tp{5})
    {
      constexpr auto s_eps = std::numeric_limits<_Tp>::epsilon();

      if (p < _Tp{0})
	return -erf_inv_recur(-p);
      else
	{
	  // Iterate experfc(x).
	  auto xprev2 = _Tp{0}, xprev = _Tp{0};
	  const auto s_max_iter = 500;
	  auto iter = 0;
	  while (++iter < s_max_iter)
	    {
	      xprev2 = xprev;
	      xprev = x;
	      x = std::log(erfc_scaled(x) / (_Tp{1} - p));
	      // If the fraction jumps < 0 just bop it back.
	      if (x < _Tp{0})
		x = -x;
	      x = std::sqrt(x);
	      if (std::abs(x - xprev) < s_eps)
		break;
	      if (std::abs(x - xprev2) < s_eps)
		break;
	    }
	  return x;
	}
    }

  /**
   * Return the inverse of the error function by series:
   * @f[
   *   erf^{-1}(x) = \sum_{k=0}^{\infty}a_k \chi^{2k+1}
   * @f]
   * where @f$ \chi = \sqrt{\pi} x / 2 @f$ and the coefficients
   * are obtained recursively:
   * @f[
   *   a_0 = 1, \hbox{    }  a_k = \frac{1}{2k+1}\sum_{j=1}^{k}
   *              \frac{a_{j-1}a_{k-j}}{j(2j-1)}
   * @f]
   */
  template<typename _Tp>
    _Tp
    erf_inv_series(_Tp p)
    {
      const auto s_sqrt_pi = _Tp{1.772453850905516027298167483341145182797L};
      const auto s_eps = std::numeric_limits<_Tp>::epsilon();

      constexpr std::size_t s_num_c = 100;
      constexpr _Tp
      s_c[s_num_c]
      {
	1.0000000000000000000000000000000000000000e+00L,
	3.3333333333333333333333333333333333333333e-01L,
	2.3333333333333333333333333333333333333333e-01L,
	2.0158730158730158730158730158730158730159e-01L,
	1.9263668430335097001763668430335097001764e-01L,
	1.9532547699214365881032547699214365881033e-01L,
	2.0593586454697565808676919788030899142010e-01L,
	2.2320975741875212774683674154573625473096e-01L,
	2.4697023314275492924730397497282028872443e-01L,
	2.7765382560322399480934889900295425728471e-01L,
	3.1614262355311719556009431370448638521126e-01L,
	3.6371758703969220002632916638901842346452e-01L,
	4.2207208084304258262391592586316200175415e-01L,
	4.9336326556393457504183570246399309984169e-01L,
	5.8029384606151398635317592916659393772989e-01L,
	6.8622339694769123804669062199054905145413e-01L,
	8.1531220555280811772727095849187818314425e-01L,
	9.7270320886455252916880107841532305635012e-01L,
	1.1647499636184417958260062794089671487378e+00L,
	1.3993010831666702410746034204406872892170e+00L,
	1.6860544545395053771174396413200857177616e+00L,
	2.0369980191940678543509421084098582912973e+00L,
	2.4669581652045463805692433503498659906947e+00L,
	2.9942820664791190057454626216924775308419e+00L,
	3.6416868900303454689338198116367649900894e+00L,
	4.4373170116450108681934433223574261955587e+00L,
	5.4160606510185387130064100587120942840981e+00L,
	6.6211901806982408845539379884204940649576e+00L,
	8.1064064311546261706653585569398214023334e+00L,
	9.9383874240100468518457098523312318792333e+00L,
	1.2199967141238320410382273076304681711955e+01L,
	1.4994101464786604505789500164531685660891e+01L,
	1.8448817910480737569763801607554579339737e+01L,
	2.2723395255097443239501906399467398597323e+01L,
	2.8016081155160551519684547291199351449624e+01L,
	3.4573733567909432191137737162607413373167e+01L,
	4.2703869211446465071921630136798549613693e+01L,
	5.2789724468807423948381760855944222262666e+01L,
	6.5309087363238988886661476478485405644500e+01L,
	8.0857851440904475210290216344379952916796e+01L,
	1.0017948355831998905863335553877479483854e+02L,
	1.2420190020854427657039294108107661830241e+02L,
	1.5408362687255315597445949256845065119843e+02L,
	1.9127159173595104118933482088635079442930e+02L,
	2.3757350384078273090471510490966501555797e+02L,
	2.9524851762450409379920541285889062493822e+02L,
	3.6712083016535879473395596350190355497002e+02L,
	4.5672204377443933850137763568856373295807e+02L,
	5.6846961707590971302659761427865879479932e+02L,
	7.0789060068197778812083087086038508771175e+02L,
	8.8190220749542727942404580881576100399193e+02L,
	1.0991637265051202368607534381089160764874e+03L,
	1.3705180086659889842849728330341658385136e+03L,
	1.7095454304103539636341488419584810137313e+03L,
	2.1332591216718800511339547506577884313893e+03L,
	2.6629776421683523690680448600112704479322e+03L,
	3.3254205938548686080054163537075449114337e+03L,
	4.1540843622567677064408565788543701928679e+03L,
	5.1909699061108619908913529600750652695848e+03L,
	6.4887530459029787592548395949069901253419e+03L,
	8.1135110231806969576401137084796791731262e+03L,
	1.0148148455685076191198091204554612819918e+04L,
	1.2696702762849418612831013617190253496734e+04L,
	1.5889755653313424153893771652199251088708e+04L,
	1.9891235834860398117901329020625469535487e+04L,
	2.4906971858173185789993829890006118493638e+04L,
	3.1195446884728408193751842686512359575872e+04L,
	3.9081324149556051730536995843174609305213e+04L,
	4.8972459241186036833087108877662536661765e+04L,
	6.1381300948238122391758767435402578263025e+04L,
	7.6951816290849023380037435887374906621876e+04L,
	9.6493370032083015879702749434102929678681e+04L,
	1.2102336029672863679530416968188136458198e+05L,
	1.5182087989802965549013218751515929183239e+05L,
	1.9049426279385044164568113689012218052878e+05L,
	2.3906611855743722117657424092133103372688e+05L,
	3.0008039496217933230379343832563534618464e+05L,
	3.7673719034285128774734288923175691261793e+05L,
	4.7306252714505085748570827450272234868736e+05L,
	5.9412217655641593142541119354284309781619e+05L,
	7.4629099296441743505320308856922156344577e+05L,
	9.3759220446465523393080671929252528124112e+05L,
	1.1781248736157831002099220455744471541354e+06L,
	1.4806024953943197978193386360494555030695e+06L,
	1.8610316947017231929650196225033400055040e+06L,
	2.3395675493065509118880041735203559358196e+06L,
	2.9415916066546777064510165052821926923954e+06L,
	3.6990707032785607548611610462302690152830e+06L,
	4.6522698884964623983695412547399497152275e+06L,
	5.8519119266215819580659637645743050898863e+06L,
	7.3619000481348488005220757750268855035060e+06L,
	9.2627511587262624888024050186222542898920e+06L,
	1.1655925260871023579609094033218903890678e+07L,
	1.4669285498580608934887412112799982558118e+07L,
	1.8463984658141537170183197023288273408861e+07L,
	2.3243151500760668256237341738654446740435e+07L,
	2.9262848204636083151389698599528628238603e+07L,
	3.6845893807083860789364230456262848641388e+07L,
	4.6399304624452616065096616113147981500063e+07L,
	5.8436299732183403140670782938460945305608e+07L,
      };

      if (p < _Tp{0})
	return -erf_inv_series(-p);
      else
	{
	  auto chi = s_sqrt_pi * p / _Tp{2};
	  auto chi2 = chi * chi;
	  auto chik = chi;
	  auto inverf = _Tp{0};
	  for (std::size_t k = 0; k < s_num_c; ++k)
	    {
	      auto term = s_c[k] * chik;
	      inverf += term;
	      if (std::abs(term) < s_eps * std::abs(inverf))
		break;
	      chik *= chi2;
	    }
	  return inverf;
	}
    }

  /**
   * Return the inverse error function.
   */
  template<typename _Tp>
    _Tp
    erf_inv(_Tp p)
    {
      const auto s_inf = emsr::infinity(p);
      if (std::isnan(p))
	return p;
      else if (std::abs(p) > _Tp(1))
	throw std::domain_error("erf_inv: Argument must have absolute value less than or equal to one.");
      else if (p == _Tp{-1})
	return -s_inf;
      else if (p == _Tp{1})
	return s_inf;
      else if (std::abs(p) > _Tp{0.95})
	return erf_inv_recur(p);
      else if (std::abs(p) > _Tp{0.75})
	return erf_inv_recur(p, erf_inv_series(p));
      else
	return erf_inv_series(p);
    }

  /**
   * Return the inverse complementary error function.
   * @todo Don't fall back on inv_erf(1-p)
   */
  template<typename _Tp>
    _Tp
    erfc_inv(_Tp q)
    {
      const auto s_inf = emsr::infinity(q);
      if (std::isnan(q))
	return q;
      else if (q < _Tp(0) || q > _Tp(2))
	throw std::domain_error("erf_inv: Argument must be within the domain of erfc: [0,2].");
      else if (q == _Tp{0})
	return +s_inf;
      else if (q == _Tp{2})
	return -s_inf;
      else if (std::abs(_Tp{1} - q) > _Tp{0.95})
	return erf_inv_recur(_Tp{1} - q);
      else if (std::abs(_Tp{1} - q) > _Tp{0.75})
	return erf_inv_recur(_Tp{1} - q, erf_inv_series(_Tp{1} - q));
      else
	return erf_inv_series(_Tp{1} - q);
    }

}

  float
  inline erf_invf(float p)
  { return detail::erf_inv<float>(p); }

  long double
  inline erf_invl(long double p)
  { return detail::erf_inv<long double>(p); }

  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    erf_inv(_Tp p)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return detail::erf_inv<type>(p);
    }

  float
  inline erfc_invf(float q)
  { return detail::erfc_inv<float>(q); }

  long double
  inline erfc_invl(long double q)
  { return detail::erfc_inv<long double>(q); }

  template<typename _Tp>
    inline emsr::fp_promote_t<_Tp>
    erfc_inv(_Tp q)
    {
      using type = emsr::fp_promote_t<_Tp>;
      return detail::erfc_inv<type>(q);
    }

/**
 * Test the inverse error function.
 */
template<typename _Tp>
  void
  test_inv_erf()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    decltype(std::cout.precision()) xw = 22;
    auto w = std::max(xw, 8 + std::cout.precision());

    const int n_max = 250;
    std::vector<_Tp> a;
    a.push_back(1);
    for (int n = 1; n < n_max; ++n)
      {
	auto atemp = _Tp{0};
	for (int k = 1; k <= n; ++k)
	  atemp += _Tp(2 * (k - 1) + 1) * a[k - 1]
		 * _Tp(2 * (n - k) + 1) * a[n - k]
		 / _Tp(k * (2 * k - 1));
	atemp /= _Tp(2 * n + 1);
	a.push_back(atemp);
      }

    std::cout << "\n\n" << std::setw(w) << " a_k" << '\n';
    for (auto aa : a)
      std::cout << ' ' << std::setw(w) << aa << '\n';

    std::vector<_Tp> c;
    c.push_back(1);
    for (int n = 1; n < n_max; ++n)
      {
	auto ctemp = _Tp{0};
	for (int k = 1; k <= n; ++k)
	  ctemp += c[k - 1] * c[n - k] / _Tp(k * (2 * k - 1));
	c.push_back(ctemp);
      }
    for (int n = 1; n < n_max; ++n)
      c[n] /= _Tp(2 * n + 1);

    std::cout << "\n\n " << std::setw(w) << "c_k" << '\n';
    for (auto cc : c)
      std::cout << ' ' << std::setw(w) << cc << '\n';

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "\"p\""
	      << ' ' << std::setw(w) << "\"inv_erf(p)\""
	      << ' ' << std::setw(w) << "\"erf(inv_erf(p))\""
	      << ' ' << std::setw(w) << "\"erf(inv_erf(p)) - p\""
	      << ' ' << std::setw(w) << "\"erf(inv_erf(p))\""
	      << ' ' << std::setw(w) << "\"erf(inv_erf(p)) - p\""
	      << ' ' << std::setw(w) << "\"erf(inv_erf(p))\""
	      << ' ' << std::setw(w) << "\"erf(inv_erf(p)) - p\""
	      << '\n';
    for (int i = -100; i <= 100; ++i)
      {
	auto p = _Tp(i * 0.01L);
	auto inverfs = detail::erf_inv_series(p);
	auto inverfr = detail::erf_inv_recur(p);
	auto inverf = erf_inv(p);
	std::cout << ' ' << std::setw(w) << p
		  << ' ' << std::setw(w) << inverfs
		  << ' ' << std::setw(w) << std::erf(inverf)
		  << ' ' << std::setw(w) << std::erf(inverf) - p
		  << ' ' << std::setw(w) << std::erf(inverfs)
		  << ' ' << std::setw(w) << std::erf(inverfs) - p
		  << ' ' << std::setw(w) << std::erf(inverfr)
		  << ' ' << std::setw(w) << std::erf(inverfr) - p
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "\"x\""
	      << ' ' << std::setw(w) << "\"erf(x)\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x))\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x)) - x\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x))\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x)) - x\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x))\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x)) - x\""
	      << '\n';
    for (int i = -200; i <= 200; ++i)
      {
	auto x = _Tp(i * 0.01L);
	auto erfx = std::erf(x);
	auto inverfs = detail::erf_inv_series(erfx);
	auto inverfr = detail::erf_inv_recur(erfx);
	auto inverf = erf_inv(erfx);
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << erfx
		  << ' ' << std::setw(w) << inverf
		  << ' ' << std::setw(w) << inverf - x
		  << ' ' << std::setw(w) << inverfs
		  << ' ' << std::setw(w) << inverfs - x
		  << ' ' << std::setw(w) << inverfr
		  << ' ' << std::setw(w) << inverfr - x
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "\"p\""
	      << ' ' << std::setw(w) << "\"inv_erfc(p)\""
	      << ' ' << std::setw(w) << "\"erfc(inv_erfc(p))\""
	      << ' ' << std::setw(w) << "\"erfc(inv_erfc(p)) - p\""
	      << '\n';
    for (int i = 200; i >= 0; --i)
      {
	auto p = _Tp(i * 0.01L);
	auto inverfc = erfc_inv(p);
	std::cout << ' ' << std::setw(w) << p
		  << ' ' << std::setw(w) << inverfc
		  << ' ' << std::setw(w) << std::erfc(inverfc)
		  << ' ' << std::setw(w) << std::erfc(inverfc) - p
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "\"x\""
	      << ' ' << std::setw(w) << "\"erf(x)\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x))\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x)) - x\""
	      << '\n';
    for (int i = -200; i <= 200; ++i)
      {
	auto x = _Tp(i * 0.01L);
	auto erfcx = std::erfc(x);
	auto inverfc = erfc_inv(erfcx);
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << erfcx
		  << ' ' << std::setw(w) << inverfc
		  << ' ' << std::setw(w) << inverfc - x
		  << '\n';
      }
  }

/**
 * Test the scaled complementary error function - experfc(x) = exp(x^2)erfc(x).
 */
template<typename _Tp>
  void
  plot_inv_erf()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    decltype(std::cout.precision()) xw = 22;
    auto w = std::max(xw, 8 + std::cout.precision());

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "p"
	      << ' ' << std::setw(w) << "inv_erf(p)"
	      << '\n';
    for (int k = -100; k <= 100; ++k)
      {
	auto p = k * _Tp{0.01L};
	auto inverf = erf_inv(p);
	std::cout << ' ' << std::setw(w) << p
		  << ' ' << std::setw(w) << inverf
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "p"
	      << ' ' << std::setw(w) << "inv_erfc(p)"
	      << '\n';
    for (int k = 200; k >= 0; --k)
      {
	auto p = k * _Tp{0.01L};
	auto inverfc = erfc_inv(p);
	std::cout << ' ' << std::setw(w) << p
		  << ' ' << std::setw(w) << inverfc
		  << '\n';
      }
  }

#ifdef STANDALONE
int
main()
{
  plot_inv_erf<double>();

  std::cout << "\n\n  float\n";
  std::cout << "  =====\n";
  test_inv_erf<float>();

  std::cout << "\n\n  double\n";
  std::cout << "  ======\n";
  test_inv_erf<double>();

  std::cout << "\n\n  long double\n";
  std::cout << "  ===========\n";
  test_inv_erf<long double>();

  //std::cout << "\n\n  __float128\n";
  //std::cout << "  ==========\n";
  //test_inv_erf<__float128>();
}
#endif
