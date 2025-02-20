/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>

#include <emsr/float128_io.h>
#include <emsr/summation.h>
#include <emsr/math_constants.h>

  constexpr unsigned long long
  s_num_harmonic_numer = 29;
  constexpr unsigned long long
  s_harmonic_numer[s_num_harmonic_numer]
  {
    1ULL,
    3ULL,
    11ULL,
    25ULL,
    137ULL,
    49ULL,
    363ULL,
    761ULL,
    7129ULL,
    7381ULL,
    83711ULL,
    86021ULL,
    1145993ULL,
    1171733ULL,
    1195757ULL,
    2436559ULL,
    42142223ULL,
    14274301ULL,
    275295799ULL,
    55835135ULL,
    18858053ULL,
    19093197ULL,
    444316699ULL,
    1347822955ULL,
    34052522467ULL,
    34395742267ULL,
    312536252003ULL,
    315404588903ULL,
    9227046511387ULL
  };
  constexpr unsigned long long
  s_harmonic_denom[s_num_harmonic_numer]
  {
    1ULL,
    2ULL,
    6ULL,
    12ULL,
    60ULL,
    20ULL,
    140ULL,
    280ULL,
    2520ULL,
    2520ULL,
    27720ULL,
    27720ULL,
    360360ULL,
    360360ULL,
    360360ULL,
    720720ULL,
    12252240ULL,
    4084080ULL,
    77597520ULL,
    15519504ULL,
    5173168ULL,
    5173168ULL,
    118982864ULL,
    356948592ULL,
    8923714800ULL,
    8923714800ULL,
    80313433200ULL,
    80313433200ULL,
    2329089562800ULL
  };

  template<typename Tp>
    Tp
    harmonic_number(unsigned int n)
    {
      if (n <= s_num_harmonic_numer)
	return Tp{s_harmonic_numer[n - 1]} / Tp{s_harmonic_denom[n - 1]};
      else
        {
	  unsigned int k = s_num_harmonic_numer - 1;
	  auto _H_k = Tp{s_harmonic_numer[k]} / Tp{s_harmonic_denom[k]};
          for (k = s_num_harmonic_numer; k <= n; ++k)
	    _H_k += Tp{1} / Tp(k);
	}
    }

  //  From sf_gamma.tcc
  template<typename Tp>
    Tp
    bernoulli_series(unsigned int n)
    {
      constexpr std::size_t num_bernoulli_numbers = 24;
      constexpr Tp
      bernoulli_numbers[num_bernoulli_numbers]
      {
		     Tp{1ULL},			-Tp{1ULL} / Tp{2ULL},
		     Tp{1ULL} /     Tp{6ULL},  Tp{0ULL},
		    -Tp{1ULL} /    Tp{30ULL},  Tp{0ULL},
		     Tp{1ULL} /    Tp{42ULL},  Tp{0ULL},
		    -Tp{1ULL} /    Tp{30ULL},  Tp{0ULL},
		     Tp{5ULL} /    Tp{66ULL},  Tp{0ULL},
		  -Tp{691ULL} /  Tp{2730ULL},  Tp{0ULL},
	             Tp{7ULL} /     Tp{6ULL},  Tp{0ULL},
	         -Tp{3617ULL} /   Tp{510ULL},  Tp{0ULL},
	         Tp{43867ULL} /   Tp{798ULL},  Tp{0ULL},
	       -Tp{174611ULL} /   Tp{330ULL},  Tp{0ULL},
	        Tp{854513ULL} /   Tp{138ULL},  Tp{0ULL}
      };

      if (n == 0)
	return Tp{1};
      else if (n == 1)
	return -Tp{1} / Tp{2};
      else if (n % 2 == 1) // Take care of the rest of the odd ones.
	return Tp{0};
      else if (n < num_bernoulli_numbers) // Return small evens that are painful for the series.
	return bernoulli_numbers[n];
      else
	{
	  Tp fact = Tp{1};
	  if ((n / 2) % 2 == 0)
	    fact *= -Tp{1};
	  for (unsigned int k = 1; k <= n; ++k)
	    fact *= k / (Tp{2} * emsr::pi_v<Tp>);
	  fact *= Tp{2};

	  Tp sum = Tp{0};
	  for (unsigned int i = 1; i < 1000; ++i)
	    {
	      Tp term = std::pow(Tp(i), -Tp(n));
	      if (term < std::numeric_limits<Tp>::epsilon())
		break;
	      sum += term;
	    }

	  return fact * sum;
	}
      return Tp{0};
    }

  /**
   *  Coefficients for Euler-Maclaurin summation of zeta functions.
   *    B_{2j} / (2j)!
   *  where B_k are the Bernoulli numbers.
   */
  constexpr size_t _Num_Euler_Maclaurin_zeta = 100;
  constexpr long double
  s_Euler_Maclaurin_zeta[_Num_Euler_Maclaurin_zeta]
  {
    1.00000000000000000000000000000000000L,
    8.33333333333333333333333333333333293e-02L,
   -1.38888888888888888888888888888888875e-03L,
    3.30687830687830687830687830687830609e-05L,
   -8.26719576719576719576719576719576597e-07L,
    2.08767569878680989792100903212014296e-08L,
   -5.28419013868749318484768220217955604e-10L,
    1.33825365306846788328269809751291227e-11L,
   -3.38968029632258286683019539124944218e-13L,
    8.58606205627784456413590545042562615e-15L,
   -2.17486869855806187304151642386591768e-16L,
    5.50900282836022951520265260890225438e-18L,
   -1.39544646858125233407076862640635480e-19L,
    3.53470703962946747169322997780379902e-21L,
   -8.95351742703754684639940801672890898e-23L,
    2.26795245233768305922449726817928506e-24L,
   -5.74479066887220244232839527972348697e-26L,
    1.45517247561486490107622443104134417e-27L,
   -3.68599494066531017606286927671534186e-29L,
    9.33673425709504466636710365024250844e-31L,
   -2.36502241570062993304902926977940878e-32L,
    5.99067176248213430064218240871649208e-34L,
   -1.51745488446829026064464819837699250e-35L,
    3.84375812545418822940606216740290214e-37L,
   -9.73635307264669102780496423687655647e-39L,
    2.46624704420068095513732412234574675e-40L,
   -6.24707674182074368796151949260113716e-42L,
    1.58240302446449142838660289637807111e-43L,
   -4.00827368594893596494573716493578672e-45L,
    1.01530758555695563022273865751378251e-46L,
   -2.57180415824187174746079460115444631e-48L,
    6.51445603523381492510893884688687796e-50L,
   -1.65013099068965245381972311645983560e-51L,
    4.17983062853947589044505904302589394e-53L,
   -1.05876346677029087587739692698831756e-54L,
    2.68187919126077066314325024787533130e-56L,
   -6.79327935110742120171687695308170512e-58L,
    1.72075776166814048850302218744398370e-59L,
   -4.35873032934889383811051833724115685e-61L,
    1.10407929036846667370868730253333297e-62L,
   -2.79666551337813450363217687544853825e-64L,
    7.08403650167947018923360554701085951e-66L,
   -1.79440740828922406419836715043711635e-67L,
    4.54528706361109610084319503460215356e-69L,
   -1.15133466319820517965514570874684656e-70L,
    2.91636477109236135051215065347539762e-72L,
   -7.38723826349733755172097357215101415e-74L,
    1.87120931176379530341655669167307802e-75L,
   -4.73982855776179939823365705874837915e-77L,
    1.20061259933545065010289119344900538e-78L,
   -3.04118724151429237818089382050570806e-80L,
    7.70341727470510626032951412805395999e-82L,
   -1.95129839090988306787181681770634078e-83L,
    4.94269656515946146653024823540482418e-85L,
   -1.25199966591718479000903037235065000e-86L,
    3.17135220176351545507490909160914066e-88L,
   -8.03312897073533444702820185587950603e-90L,
    2.03481533916614656707738578184657923e-91L,
   -5.15424746644747384952223952829841139e-93L,
    1.30558613521494672211652811590429162e-94L,
   -3.30708831417509124211476473245870569e-96L,
    8.37695256004909128671576487001515732e-98L,
   -2.12190687174971376532740302523869392e-99L,
    5.37485289561228024546639030950850747e-101L,
   -1.36146614321720693646766878619458764e-102L,
    3.44863402799339902711847245019887448e-104L,
   -8.73549204163835504185857147126132564e-106L,
    2.21272598339254970369646016351296266e-107L,
   -5.60490039283722415865004549966830368e-109L,
    1.41973785499917876418113219390120402e-110L,
   -3.59623799825876265563506189711913868e-112L,
    9.10937726607823184392152343960921788e-114L,
   -2.30743221710912328319788632553096580e-115L,
    5.84479408529900198086678715736051924e-117L,
   -1.48050363717057449175295580276805879e-118L,
    3.75015952262271968009868219553161640e-120L,
   -9.49926504199295827433414434212738749e-122L,
    2.40619194446751986855735721071770656e-123L,
   -6.09495539710268473564828955326838575e-125L,
    1.54387023770424714174247152273407761e-126L,
   -3.91066899685929231015496868644244724e-128L,
    9.90584028987942974230798861783873238e-130L,
   -2.50917865785635531226030684021698888e-131L,
    6.35582378960245978851062983907498289e-133L,
   -1.60994897346162503541124270943910111e-134L,
    4.07804838987246223591940353095558419e-136L,
   -1.03298172453151901860990367375763718e-137L,
    2.61657327526092016567608044513104821e-139L,
   -6.62785753340862278242754283133858879e-141L,
    1.67885592574436730038021560419761003e-142L,
   -4.25259173903430060933988427590406140e-144L,
    1.07719407136645728953530432238003752e-145L,
   -2.72856445808395795763235897809610861e-147L,
    6.91153451343701367068428806267885056e-149L,
   -1.75071214421576824424282262560183841e-150L,
    4.43460566672391957275357495883452114e-152L,
   -1.12329873784871496672434722468912147e-153L,
    2.84534894256939708475369589042384123e-155L,
   -7.20735306841513677412535382503699916e-157L,
    1.82564385955014175253212078464905862e-158L
  };


  template<typename Tp>
    Tp
    hurwitz_zeta_euler_maclaurin(Tp s, Tp a)
    {
      constexpr auto s_eps = std::numeric_limits<Tp>::epsilon();
      constexpr auto s_N = 10 + std::numeric_limits<Tp>::digits10 / 2;
      constexpr auto s_jmax = 99;

      const auto pmax  = std::pow(Tp(s_N) + a, -s);
      const auto denom = (s_N + a) * (s_N + a);
      auto ans = pmax * ((s_N + a) / (s - Tp{1}) + Tp{0.5L});
      for(auto k = 0; k < s_N; ++k)
        ans += std::pow(k + a, -s);

      auto fact = pmax * s / (s_N + a);
      auto delta_prev = std::numeric_limits<Tp>::max();
      for(auto j = 0; j <= s_jmax; ++j)
        {
	  auto delta = s_Euler_Maclaurin_zeta[j + 1] * fact;
	  if (std::abs(delta) > delta_prev)
	    break;
	  delta_prev = std::abs(delta);
	  ans += delta;
	  if(std::abs(delta / ans) < Tp{0.5L} * s_eps)
	    break;
	  fact *= (s + Tp(2 * j + 1)) * (s + Tp(2 * j + 2))
		  / denom;
        }

      return ans;
    }

  template<typename Tp>
    Tp
    riemann_zeta_euler_maclaurin(Tp s)
    {
      constexpr auto s_eps = std::numeric_limits<Tp>::epsilon();
      constexpr auto s_N = 10 + std::numeric_limits<Tp>::digits10 / 2;
      constexpr auto s_jmax = 99;

      const auto pmax  = std::pow(Tp(s_N) + Tp{1}, -s);
      const auto denom = (s_N + Tp{1}) * (s_N + Tp{1});
      auto ans = pmax * ((s_N + Tp{1}) / (s - Tp{1}) + Tp{0.5L});
      for(auto k = 0; k < s_N; ++k)
        ans += std::pow(k + Tp{1}, -s);

      auto fact = pmax * s / (s_N + Tp{1});
      auto delta_prev = std::numeric_limits<Tp>::max();
      for(auto j = 0; j <= s_jmax; ++j)
        {
	  auto delta = s_Euler_Maclaurin_zeta[j + 1] * fact;
	  if (std::abs(delta) > delta_prev)
	    break;
	  delta_prev = std::abs(delta);
	  ans += delta;
	  if(std::abs(delta / ans) < Tp{0.5L} * s_eps)
	    break;
	  fact *= (s + Tp(2 * j + 1)) * (s + Tp(2 * j + 2))
		  / denom;
        }

      return ans;
    }

  template<typename Tp>
    Tp
    hurwitz_zeta_glob(Tp s, Tp a)
    {
      const auto s_eps = std::numeric_limits<Tp>::epsilon();
      // Max before overflow?
      const auto s_max = std::numeric_limits<Tp>::max();
      const auto s_inf = std::numeric_limits<Tp>::infinity();

      //std::cout.precision(std::numeric_limits<Tp>::max_digits10);

      if (s == +Tp{0})
	return s_inf;

      std::vector<Tp> apow;

      constexpr unsigned int s_maxit = 1000;
      // Zeroth order contribution already calculated.
      const auto a1ms = std::pow(a, Tp{1} - s);
      apow.push_back(a1ms);
      auto zeta = Tp(a1ms);
#ifdef DEBUG_SERIES
      std::cout << "    n=" << 0 << " term=" << zeta << " zeta=" << zeta << '\n';
#endif
      for (unsigned int n = 1; n < s_maxit; ++n)
	{
	  bool punt = false;
	  // Again, the zeroth order contribution already calculated.
	  auto termp = Tp(a1ms);
	  auto termm = Tp{0};
	  auto binom = Tp{1};
#ifdef DEBUG_SERIES
	  std::cout << "        n=" << setw(4) << n << " k=" << setw(4) << 0
		    << " binom=" << setw(20) << binom
		    << " termp=" << setw(20) << termp
		    << " termm=" << setw(20) << termm << '\n';
#endif
	  apow.push_back(std::pow(Tp(n) + a, Tp{1} - s));
	  for (unsigned int k = 1; k <= n; ++k)
	    {
	      binom *= Tp(n - k + 1) / Tp(k);
	      if (std::abs(binom) > s_max)
		{
		  punt = true;
		  break;
		}
	      (k % 2 == 0 ? termp : termm) += binom * apow[k];
#ifdef DEBUG_SERIES
	      std::cout << "        n=" << setw(4) << n << " k=" << setw(4) << k
			<< " binom=" << setw(20) << binom
			<< " termp=" << setw(20) << termp
			<< " termm=" << setw(20) << termm << '\n';
#endif
	    }
	  if (punt)
	    break;
	  auto term = (termp - termm) / (n + 1);
	  zeta += term;
	  if (std::abs(term) < s_eps * std::abs(zeta))
	    break;
#ifdef DEBUG_SERIES
	  std::cout << "    n=" << n << " term=" << term << " zeta=" << zeta << '\n';
#endif
	}

      zeta /= s - Tp{1};

      return zeta;
    }

  //  I bet you could analytically continue this with gamma(Tp(m + 1))
  template<typename Tp>
    Tp
    polygamma(unsigned int m, Tp z)
    {
      auto sign = (m % 2 == 0 ? -1 : +1);
      Tp factorial{1};
      for (unsigned int k = 1; k <= m; ++k)
	factorial *= k;
      return sign * factorial * hurwitz_zeta_euler_maclaurin(Tp(m + 1), z);
    }


  template<typename Tp>
    Tp
    riemann_zeta_m_1_basic_sum(Tp s)
    {
      constexpr auto s_eps = emsr::epsilon<Tp>();
      // A user shouldn't get to this.
      if (s < Tp{1})
	throw std::domain_error("Bad argument in zeta sum.");

      int k_max = std::min(1000000, int(std::pow(Tp{1} / s_eps, Tp{1} / s)));
#ifdef DEBUG_SERIES
      std::cerr << "s = " << s << "  k_max = " << k_max << '\n';
#endif
      auto zeta_m_1 = Tp{0};
      for (int k = k_max; k >= 2; --k)
	{
	  auto term = std::pow(Tp(k), -s);
	  zeta_m_1 += term;
	  if (term < s_eps * zeta_m_1)
	    break;
	}

      return zeta_m_1;
    }


  template<typename Tp>
    Tp
    riemann_zeta_m_1_vanwg_sum(Tp s)
    {
      constexpr auto s_eps = emsr::epsilon<Tp>();
      // A user shouldn't get to this.
      if (s < Tp{1})
	throw std::domain_error("Bad argument in zeta sum.");

      int k_max = std::min(1000000, int(std::pow(Tp{1} / s_eps, Tp{1} / s)));
#ifdef DEBUG_SERIES
      std::cerr << "s = " << s << "  k_max = " << k_max << '\n';
#endif
      auto zeta_m_1 = Tp{0};
      for (int k = k_max; k >= 2; --k)
	{
	  auto term = std::pow(Tp(k), -s);
	  zeta_m_1 += term;
	  if (term < s_eps * zeta_m_1)
	    break;
	}

      return zeta_m_1;
    }


  template<typename Tp>
    Tp
    riemann_zeta_m_1_kahan_sum(Tp s)
    {
      constexpr auto s_eps = emsr::epsilon<Tp>();
      // A user shouldn't get to this.
      if (s < Tp{1})
	throw std::domain_error("Bad argument in zeta sum.");

      emsr::KahanSum<Tp> zeta_m_1;
      int k_max = std::min(1000000, int(std::pow(Tp{1} / s_eps, Tp{1} / s)));
#ifdef DEBUG_SERIES
      std::cerr << "s = " << s << "  k_max = " << k_max << '\n';
#endif
      for (int k = k_max; k >= 2; --k)
      //int k_max = 10000;
      //for (int k = k_max; k >= 2; --k)
	{
	  auto term = std::pow(Tp(k), -s);
	  zeta_m_1 += term;
	  if (term < s_eps * zeta_m_1())
	    break;
	}

      return zeta_m_1();
    }


template<typename Tp>
  void
  test_hurwitz_zeta_new()
  {
    std::cout.precision(std::numeric_limits<Tp>::max_digits10);
    auto width = std::numeric_limits<Tp>::max_digits10 + 6;

    // Build the harmonic numbers H_n.
    std::cout << "\nBuild the H_n numbers\n";
    auto hn = Tp{0};
    for (auto i = 1; i < 100; ++i)
      {
	std::cout << (hn += Tp{1} / i) << '\n';
      }

    // Build the B_{2j} numbers.
    std::cout << "\nBuild the B_{2j} numbers\n";
    for (auto i = 1; i < 100; ++i)
      {
	std::cout << bernoulli_series<Tp>(2 * i) << '\n';
      }

    // Build the B_{2j}/(2j)! numbers.
    std::cout << "\nBuild the B_{2j}/(2j)! numbers\n";
    Tp fact{1};
    for (auto i = 1; i < 100; ++i)
      {
	fact /= (2 * i - 1) * (2 * i);
	std::cout << fact * bernoulli_series<Tp>(2 * i) << '\n';
      }

    // Test zeta - 1 function with both simple and Kahan summation.
    std::cout << "\nTest zeta - 1 function with both simple and Kahan summation\n";
    for (auto is = 10; is < 100; ++is)
      {
	Tp s = 0.1L * is;
	auto zetam1s = riemann_zeta_m_1_basic_sum(s);
	auto zetam1k = riemann_zeta_m_1_kahan_sum(s);
	std::cout << ' ' << std::setw(width) << s
		  << ' ' << std::setw(width) << zetam1s
		  << ' ' << std::setw(width) << zetam1k
		  << ' ' << std::setw(width) << zetam1k - zetam1s
		  << '\n';
      }

    // Test zeta - 1 function with both simple and Kahan summation.
    std::cout << "\nTest zeta - 1 function with both simple and Kahan summation\n";
    for (auto is = 1; is <= 100; ++is)
      {
	Tp s = 1.0L * is;
	auto zetam1s = riemann_zeta_m_1_basic_sum(s);
	auto zetam1k = riemann_zeta_m_1_kahan_sum(s);
	std::cout << ' ' << std::setw(width) << s
		  << ' ' << std::setw(width) << zetam1s
		  << ' ' << std::setw(width) << zetam1k
		  << ' ' << std::setw(width) << zetam1k - zetam1s
		  << '\n';
      }

    // Test a Bernoulli thing for the regular zeta function.
    std::cout << "\nBernoulli sum regular zeta function\n";
    for (auto is = -9; is < 100; ++is)
      {
	Tp s = 0.1L * is;
	auto hzeta = hurwitz_zeta_euler_maclaurin(s, Tp{1});
	auto rzeta = riemann_zeta_euler_maclaurin(s);
	std::cout << ' ' << std::setw(width) << s
		  << ' ' << std::setw(width) << hzeta
		  << ' ' << std::setw(width) << rzeta
		  << ' ' << std::setw(width) << hzeta - rzeta
		  << '\n';
      }

    // Test a Hurwitz zeta function.
    std::cout << "\nHurwitz zeta function\n";
    std::cout << "\n a = " << Tp{4} << '\n';//'\n';
    std::cout << ' ' << std::setw(width) << Tp{10}
	      << ' ' << std::setw(width) << hurwitz_zeta_euler_maclaurin(Tp{10}, Tp{4})
	      << ' ' << std::setw(width) << hurwitz_zeta_glob(Tp{10}, Tp{4})
	      << '\n';

    std::cout << "\nHurwitz zeta function\n";
    for (auto ia = 1; ia < 100; ++ia)
      {
	Tp a = 0.1L * ia;
	std::cout << "\n a = " << a << '\n';//'\n';
	for (auto is = 0; is < 100; ++is)
	  {
	    Tp s = 0.1L * is;
	    if (s == 1)
	      continue;
	    std::cout << ' ' << std::setw(width) << s
		      << ' ' << std::setw(width) << hurwitz_zeta_euler_maclaurin(s, a)
		      << ' ' << std::setw(width) << hurwitz_zeta_glob(s, a)
		      << '\n';
	  }
      }
  }

int
main()
{
  test_hurwitz_zeta_new<long double>();
}

