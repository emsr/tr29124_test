/**
 *
 */

#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

#include <emsr/float128_io.h>
#include <emsr/math_constants.h>

  //  From sf_gamma.tcc
  template<typename Tp>
    Tp
    bernoulli_series(unsigned int n)
    {
      constexpr std::size_t num_bernoulli_numbers = 24;
      constexpr Tp
      bernoulli_numbers[num_bernoulli_numbers]
      {
	 Tp{1UL},	                 -Tp{1UL} / Tp{2UL},
	 Tp{1UL} / Tp{6UL},             Tp{0UL},
	-Tp{1UL} / Tp{30UL},            Tp{0UL},
	 Tp{1UL} / Tp{42UL},            Tp{0UL},
	-Tp{1UL} / Tp{30UL},            Tp{0UL},
	 Tp{5UL} / Tp{66UL},            Tp{0UL},
	-Tp{691UL} / Tp{2730UL},        Tp{0UL},
	 Tp{7UL} / Tp{6UL},             Tp{0UL},
	-Tp{3617UL} / Tp{510UL},        Tp{0UL},
	 Tp{43867UL} / Tp{798UL},       Tp{0UL},
	-Tp{174611UL} / Tp{330UL},      Tp{0UL},
	 Tp{854513UL} / Tp{138UL},      Tp{0UL}
      };

      if (n == 0)
	return Tp{1};
      else if (n == 1)
	return -Tp{1} / Tp{2};
      else if (n % 2 == 1) // Take care of the rest of the odd ones.
	return Tp{0};
      else if (n < num_bernoulli_numbers)
	// Return small evens that are painful for the series.
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
    }

  template<typename Tp>
    Tp
    hurwitz_zeta_euler_maclaurin(Tp s, Tp a)
    {
      constexpr auto s_eps = std::numeric_limits<Tp>::epsilon();
      constexpr int s_N = 10 + std::numeric_limits<Tp>::digits10 / 2;
      constexpr int s_jmax = 99;

      //  Coefficients for Euler-Maclaurin summation:
      //  B_{2j}/(2j)!
      static constexpr Tp
      s_hzeta_c[100]
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

      const auto pmax  = std::pow(Tp(s_N) + a, -s);
      const auto denom = (s_N + a) * (s_N + a);
      auto ans = pmax * ((s_N + a) / (s - Tp{1}) + Tp{0.5L});
      for(auto k = 0; k < s_N; ++k)
        ans += std::pow(k + a, -s);

      auto fact = pmax * s / (s_N + a);
      auto delta_prev = std::numeric_limits<Tp>::max();
      for(auto j = 0; j < s_jmax; ++j)
        {
	  auto delta = s_hzeta_c[j + 1] * fact;
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

      constexpr unsigned int s_maxit = 10000;
       //  Zeroth order contribution already calculated.
      auto zeta = Tp{0.5L};
      for (unsigned int n = 1; n < s_maxit; ++n)
	{
	  bool punt = false;
	  auto term = Tp{1}; // Again, the zeroth order.
	  auto binom = Tp{1};
	  for (unsigned int k = 1; k <= n; ++k)
	    {
	      binom *= -Tp(n - k + 1) / Tp(k);
	      if (std::abs(binom) > s_max)
	      {
		punt = true;
		break;
	      }
	      term += binom * std::pow(Tp(k + a), Tp{1} - s);
	      //std::cout << "        " << k << ' ' << binom << ' ' << term << '\n';
	    }
	  if (punt)
	    break;
	  term /= Tp(n + 1);
	  zeta += term;
	  if (std::abs(term / zeta) < s_eps)
	    break;
	  //std::cout << "    " << n << ' ' << term << ' ' << zeta << '\n';
	}

      zeta /= s - Tp{1};

      return zeta;
    }

template<typename Tp>
  void
  test_hurwitz_zeta()
  {
    std::cout.precision(std::numeric_limits<Tp>::max_digits10);
    auto width = std::numeric_limits<Tp>::max_digits10 + 8;

    Tp fact{1};
    for (auto i = 1; i < 100; ++i)
      {
	fact /= (2 * i - 1) * (2 * i);
	std::cout << fact * bernoulli_series<Tp>(2 * i) << '\n';
      }

    std::cout << '\n';
    for (auto ia = 0; ia < 100; ++ia)
      {
	Tp a = 0.1L * ia;
	std::cout << "a = " << a << '\n';
	for (auto is = 0; is < 100; ++is)
	  {
	    Tp s = 0.1L * is;
	    if (s == 1)
	      continue;
	    std::cout << ' ' << std::setw(4) << s
		      << ' ' << std::setw(width) << hurwitz_zeta_euler_maclaurin(s, a)
		      //<< ' ' << std::setw(width) << hurwitz_zeta_glob(s, a)
		      << '\n';
	  }
      }
  }

int
main()
{
  test_hurwitz_zeta<long double>();
}

