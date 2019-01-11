/*
$HOME/bin/bin/g++ -std=gnu++2a -DSTANDALONE -g -Wall -Wextra -Wno-psabi -I. -o test_inv_erf test_inv_erf.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_inv_erf > test_inv_erf.txt

$HOME/bin/bin/g++ -std=gnu++2a -DSTANDALONE -g -Wall -Wextra -Wno-psabi -I. -o test_inv_erf test_inv_erf.cpp -lquadmath -lmpfr
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_inv_erf > test_inv_erf.txt
*/

// AAOF pp. 408-409.

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>
#include <bits/float128_io.h>

  template<typename _Tp>
    _Tp
    __erfc_scaled(_Tp __x)
    { return std::exp(__x * __x) * std::erfc(__x); }

  /**
   * Return the inverse error function by recursion:
   * @f[
   *
   * @f]
   *
   * @param[in] __p The argument between -1 and 1
   * @param[in] __x The initial x-value guess.
   */
  template<typename _Tp>
    _Tp
    __erf_inv_recur(_Tp __p, _Tp __x = _Tp{5})
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      if (__p < _Tp{0})
	return -__erf_inv_recur(-__p);
      else
	{
	  // Iterate experfc(x).
	  auto __xprev2 = _Tp{0}, __xprev = _Tp{0};
	  const auto _S_max_iter = 500;
	  auto __iter = 0;
	  while (++__iter < _S_max_iter)
	    {
	      __xprev2 = __xprev;
	      __xprev = __x;
	      __x = std::log(__erfc_scaled(__x) / (_Tp{1} - __p));
	      // If the fraction jumps < 0 just bop it back.
	      if (__x < _Tp{0})
		__x = -__x;
	      __x = std::sqrt(__x);
	      if (std::abs(__x - __xprev) < _S_eps)
		break;
	      if (std::abs(__x - __xprev2) < _S_eps)
		break;
	    }
	  return __x;
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
    __erf_inv_series(_Tp __p)
    {
      const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797Q};
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();

      constexpr std::size_t _S_num_c = 100;
      constexpr _Tp
      _S_c[_S_num_c]
      {
	1.0000000000000000000000000000000000000000e+00Q,
	3.3333333333333333333333333333333333333333e-01Q,
	2.3333333333333333333333333333333333333333e-01Q,
	2.0158730158730158730158730158730158730159e-01Q,
	1.9263668430335097001763668430335097001764e-01Q,
	1.9532547699214365881032547699214365881033e-01Q,
	2.0593586454697565808676919788030899142010e-01Q,
	2.2320975741875212774683674154573625473096e-01Q,
	2.4697023314275492924730397497282028872443e-01Q,
	2.7765382560322399480934889900295425728471e-01Q,
	3.1614262355311719556009431370448638521126e-01Q,
	3.6371758703969220002632916638901842346452e-01Q,
	4.2207208084304258262391592586316200175415e-01Q,
	4.9336326556393457504183570246399309984169e-01Q,
	5.8029384606151398635317592916659393772989e-01Q,
	6.8622339694769123804669062199054905145413e-01Q,
	8.1531220555280811772727095849187818314425e-01Q,
	9.7270320886455252916880107841532305635012e-01Q,
	1.1647499636184417958260062794089671487378e+00Q,
	1.3993010831666702410746034204406872892170e+00Q,
	1.6860544545395053771174396413200857177616e+00Q,
	2.0369980191940678543509421084098582912973e+00Q,
	2.4669581652045463805692433503498659906947e+00Q,
	2.9942820664791190057454626216924775308419e+00Q,
	3.6416868900303454689338198116367649900894e+00Q,
	4.4373170116450108681934433223574261955587e+00Q,
	5.4160606510185387130064100587120942840981e+00Q,
	6.6211901806982408845539379884204940649576e+00Q,
	8.1064064311546261706653585569398214023334e+00Q,
	9.9383874240100468518457098523312318792333e+00Q,
	1.2199967141238320410382273076304681711955e+01Q,
	1.4994101464786604505789500164531685660891e+01Q,
	1.8448817910480737569763801607554579339737e+01Q,
	2.2723395255097443239501906399467398597323e+01Q,
	2.8016081155160551519684547291199351449624e+01Q,
	3.4573733567909432191137737162607413373167e+01Q,
	4.2703869211446465071921630136798549613693e+01Q,
	5.2789724468807423948381760855944222262666e+01Q,
	6.5309087363238988886661476478485405644500e+01Q,
	8.0857851440904475210290216344379952916796e+01Q,
	1.0017948355831998905863335553877479483854e+02Q,
	1.2420190020854427657039294108107661830241e+02Q,
	1.5408362687255315597445949256845065119843e+02Q,
	1.9127159173595104118933482088635079442930e+02Q,
	2.3757350384078273090471510490966501555797e+02Q,
	2.9524851762450409379920541285889062493822e+02Q,
	3.6712083016535879473395596350190355497002e+02Q,
	4.5672204377443933850137763568856373295807e+02Q,
	5.6846961707590971302659761427865879479932e+02Q,
	7.0789060068197778812083087086038508771175e+02Q,
	8.8190220749542727942404580881576100399193e+02Q,
	1.0991637265051202368607534381089160764874e+03Q,
	1.3705180086659889842849728330341658385136e+03Q,
	1.7095454304103539636341488419584810137313e+03Q,
	2.1332591216718800511339547506577884313893e+03Q,
	2.6629776421683523690680448600112704479322e+03Q,
	3.3254205938548686080054163537075449114337e+03Q,
	4.1540843622567677064408565788543701928679e+03Q,
	5.1909699061108619908913529600750652695848e+03Q,
	6.4887530459029787592548395949069901253419e+03Q,
	8.1135110231806969576401137084796791731262e+03Q,
	1.0148148455685076191198091204554612819918e+04Q,
	1.2696702762849418612831013617190253496734e+04Q,
	1.5889755653313424153893771652199251088708e+04Q,
	1.9891235834860398117901329020625469535487e+04Q,
	2.4906971858173185789993829890006118493638e+04Q,
	3.1195446884728408193751842686512359575872e+04Q,
	3.9081324149556051730536995843174609305213e+04Q,
	4.8972459241186036833087108877662536661765e+04Q,
	6.1381300948238122391758767435402578263025e+04Q,
	7.6951816290849023380037435887374906621876e+04Q,
	9.6493370032083015879702749434102929678681e+04Q,
	1.2102336029672863679530416968188136458198e+05Q,
	1.5182087989802965549013218751515929183239e+05Q,
	1.9049426279385044164568113689012218052878e+05Q,
	2.3906611855743722117657424092133103372688e+05Q,
	3.0008039496217933230379343832563534618464e+05Q,
	3.7673719034285128774734288923175691261793e+05Q,
	4.7306252714505085748570827450272234868736e+05Q,
	5.9412217655641593142541119354284309781619e+05Q,
	7.4629099296441743505320308856922156344577e+05Q,
	9.3759220446465523393080671929252528124112e+05Q,
	1.1781248736157831002099220455744471541354e+06Q,
	1.4806024953943197978193386360494555030695e+06Q,
	1.8610316947017231929650196225033400055040e+06Q,
	2.3395675493065509118880041735203559358196e+06Q,
	2.9415916066546777064510165052821926923954e+06Q,
	3.6990707032785607548611610462302690152830e+06Q,
	4.6522698884964623983695412547399497152275e+06Q,
	5.8519119266215819580659637645743050898863e+06Q,
	7.3619000481348488005220757750268855035060e+06Q,
	9.2627511587262624888024050186222542898920e+06Q,
	1.1655925260871023579609094033218903890678e+07Q,
	1.4669285498580608934887412112799982558118e+07Q,
	1.8463984658141537170183197023288273408861e+07Q,
	2.3243151500760668256237341738654446740435e+07Q,
	2.9262848204636083151389698599528628238603e+07Q,
	3.6845893807083860789364230456262848641388e+07Q,
	4.6399304624452616065096616113147981500063e+07Q,
	5.8436299732183403140670782938460945305608e+07Q,
      };

      if (__p < _Tp{0})
	return -__erf_inv_series(-__p);
      else
	{
	  auto __chi = _S_sqrt_pi * __p / _Tp{2};
	  auto __chi2 = __chi * __chi;
	  auto __chik = __chi;
	  auto __inverf = _Tp{0};
	  for (std::size_t __k = 0; __k < _S_num_c; ++__k)
	    {
	      auto __term = _S_c[__k] * __chik;
	      __inverf += __term;
	      if (std::abs(__term) < _S_eps * std::abs(__inverf))
		break;
	      __chik *= __chi2;
	    }
	  return __inverf;
	}
    }

  /**
   * Return the inverse error function.
   */
  template<typename _Tp>
    _Tp
    __erf_inv(_Tp __p)
    {
      const auto _S_inf = __gnu_cxx::__infinity(__p);
      if (std::isnan(__p))
	return __p;
      else if (std::abs(__p) > _Tp(1))
	std::__throw_domain_error("__erf_inv: Argument must have absolute value"
				  " less than or equal to one.");
      else if (__p == _Tp{-1})
	return -_S_inf;
      else if (__p == _Tp{1})
	return _S_inf;
      else if (std::abs(__p) > _Tp{0.95})
	return __erf_inv_recur(__p);
      else if (std::abs(__p) > _Tp{0.75})
	return __erf_inv_recur(__p, __erf_inv_series(__p));
      else
	return __erf_inv_series(__p);
    }

  /**
   * Return the inverse complementary error function.
   * @todo Don't fall back on inv_erf(1-p)
   */
  template<typename _Tp>
    _Tp
    __erfc_inv(_Tp __q)
    {
      const auto _S_inf = __gnu_cxx::__infinity(__q);
      if (std::isnan(__q))
	return __q;
      else if (__q < _Tp(0) || __q > _Tp(2))
	std::__throw_domain_error("__erf_inv: Argument must be within"
				  " the domain of erfc: [0,2].");
      else if (__q == _Tp{0})
	return +_S_inf;
      else if (__q == _Tp{2})
	return -_S_inf;
      else if (std::abs(_Tp{1} - __q) > _Tp{0.95})
	return __erf_inv_recur(_Tp{1} - __q);
      else if (std::abs(_Tp{1} - __q) > _Tp{0.75})
	return __erf_inv_recur(_Tp{1} - __q, __erf_inv_series(_Tp{1} - __q));
      else
	return __erf_inv_series(_Tp{1} - __q);
    }

  float
  inline erf_invf(float __p)
  { return __erf_inv<float>(__p); }

  long double
  inline erf_invl(long double __p)
  { return __erf_inv<long double>(__p); }

  template<typename _Tp>
    inline __gnu_cxx::fp_promote_t<_Tp>
    erf_inv(_Tp __p)
    {
      using __type = __gnu_cxx::fp_promote_t<_Tp>;
      return __erf_inv<__type>(__p);
    }

  float
  inline erfc_invf(float __q)
  { return __erfc_inv<float>(__q); }

  long double
  inline erfc_invl(long double __q)
  { return __erfc_inv<long double>(__q); }

  template<typename _Tp>
    inline __gnu_cxx::fp_promote_t<_Tp>
    erfc_inv(_Tp __q)
    {
      using __type = __gnu_cxx::fp_promote_t<_Tp>;
      return __erfc_inv<__type>(__q);
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
    for (int __n = 1; __n < n_max; ++__n)
      {
	auto __atemp = _Tp{0};
	for (int __k = 1; __k <= __n; ++__k)
	  __atemp += _Tp(2 * (__k - 1) + 1) * a[__k - 1]
		 * _Tp(2 * (__n - __k) + 1) * a[__n - __k]
		 / _Tp(__k * (2 * __k - 1));
	__atemp /= _Tp(2 * __n + 1);
	a.push_back(__atemp);
      }

    std::cout << "\n\n" << std::setw(w) << " a_k" << '\n';
    for (auto __aa : a)
      std::cout << ' ' << std::setw(w) << __aa << '\n';

    std::vector<_Tp> c;
    c.push_back(1);
    for (int __n = 1; __n < n_max; ++__n)
      {
	auto __ctemp = _Tp{0};
	for (int __k = 1; __k <= __n; ++__k)
	  __ctemp += c[__k - 1] * c[__n - __k] / _Tp(__k * (2 * __k - 1));
	c.push_back(__ctemp);
      }
    for (int __n = 1; __n < n_max; ++__n)
      c[__n] /= _Tp(2 * __n + 1);

    std::cout << "\n\n " << std::setw(w) << "c_k" << '\n';
    for (auto __cc : c)
      std::cout << ' ' << std::setw(w) << __cc << '\n';

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
    for (int __i = -100; __i <= 100; ++__i)
      {
	auto __p = _Tp(__i * 0.01Q);
	auto __inverfs = __erf_inv_series(__p);
	auto __inverfr = __erf_inv_recur(__p);
	auto __inverf = __erf_inv(__p);
	std::cout << ' ' << std::setw(w) << __p
		  << ' ' << std::setw(w) << __inverfs
		  << ' ' << std::setw(w) << std::erf(__inverf)
		  << ' ' << std::setw(w) << std::erf(__inverf) - __p
		  << ' ' << std::setw(w) << std::erf(__inverfs)
		  << ' ' << std::setw(w) << std::erf(__inverfs) - __p
		  << ' ' << std::setw(w) << std::erf(__inverfr)
		  << ' ' << std::setw(w) << std::erf(__inverfr) - __p
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
    for (int __i = -200; __i <= 200; ++__i)
      {
	auto __x = _Tp(__i * 0.01Q);
	auto __erfx = std::erf(__x);
	auto __inverfs = __erf_inv_series(__erfx);
	auto __inverfr = __erf_inv_recur(__erfx);
	auto __inverf = __erf_inv(__erfx);
	std::cout << ' ' << std::setw(w) << __x
		  << ' ' << std::setw(w) << __erfx
		  << ' ' << std::setw(w) << __inverf
		  << ' ' << std::setw(w) << __inverf - __x
		  << ' ' << std::setw(w) << __inverfs
		  << ' ' << std::setw(w) << __inverfs - __x
		  << ' ' << std::setw(w) << __inverfr
		  << ' ' << std::setw(w) << __inverfr - __x
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "\"p\""
	      << ' ' << std::setw(w) << "\"inv_erfc(p)\""
	      << ' ' << std::setw(w) << "\"erfc(inv_erfc(p))\""
	      << ' ' << std::setw(w) << "\"erfc(inv_erfc(p)) - p\""
	      << '\n';
    for (int __i = 200; __i >= 0; --__i)
      {
	auto __p = _Tp(__i * 0.01Q);
	auto __inverfc = __erfc_inv(__p);
	std::cout << ' ' << std::setw(w) << __p
		  << ' ' << std::setw(w) << __inverfc
		  << ' ' << std::setw(w) << std::erfc(__inverfc)
		  << ' ' << std::setw(w) << std::erfc(__inverfc) - __p
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "\"x\""
	      << ' ' << std::setw(w) << "\"erf(x)\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x))\""
	      << ' ' << std::setw(w) << "\"inv_erf(erf(x)) - x\""
	      << '\n';
    for (int __i = -200; __i <= 200; ++__i)
      {
	auto __x = _Tp(__i * 0.01Q);
	auto __erfcx = std::erfc(__x);
	auto __inverfc = __erfc_inv(__erfcx);
	std::cout << ' ' << std::setw(w) << __x
		  << ' ' << std::setw(w) << __erfcx
		  << ' ' << std::setw(w) << __inverfc
		  << ' ' << std::setw(w) << __inverfc - __x
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
    for (int __k = -100; __k <= 100; ++__k)
      {
	auto __p = __k * _Tp{0.01Q};
	auto __inverf = __erf_inv(__p);
	std::cout << ' ' << std::setw(w) << __p
		  << ' ' << std::setw(w) << __inverf
		  << '\n';
      }

    std::cout << "\n\n"
	      << ' ' << std::setw(w) << "p"
	      << ' ' << std::setw(w) << "inv_erfc(p)"
	      << '\n';
    for (int __k = 200; __k >= 0; --__k)
      {
	auto __p = __k * _Tp{0.01Q};
	auto __inverfc = __erfc_inv(__p);
	std::cout << ' ' << std::setw(w) << __p
		  << ' ' << std::setw(w) << __inverfc
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

  std::cout << "\n\n  __float128\n";
  std::cout << "  ==========\n";
  test_inv_erf<__float128>();
}
#endif
