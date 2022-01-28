/**
 *
 */

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <emsr/float128_io.h>
#include <emsr/numeric_limits_float128.h>
#include <emsr/math_constants.h>

  //  Checbyshev coefficient matrix.
  constexpr int
  cheby(unsigned int n, unsigned int k)
  {
    if (k > n)
      return 0;
    else if (n == 1)
      return 1;
    else if (n == 2)
      {
	if (k == 1)
	  return 0;
	else if (k == 2)
	  return 1;
      }
    else
      {
	if (k == 1)
	  return -cheby(n - 2, 1);
	else if (k == n)
	  return 2 * cheby(n - 1, k - 1);
	else
	  return 2 * cheby(n - 1, k - 1) - cheby(n - 2, k);
      }
    return -1;
  }

  template<typename Tp>
    Tp
    p(int k, Tp g)
    {
      const auto s_pi  = emsr::pi_v<Tp>;
      auto fact = std::sqrt(Tp{2} / s_pi);
      auto sum = cheby(2 * k + 1, 1) * fact
		 * std::exp(Tp(g + 0.5Q))
		 / std::sqrt(Tp(g + 0.5Q));
      for (int a = 1; a <= k; ++a)
	{
	  fact *= Tp(2 * a - 1) / 2;
	  sum += cheby(2 * k + 1, 2 * a + 1) * fact
		 * std::pow(Tp(a + g + 0.5Q), -Tp(a + 0.5Q))
		 * std::exp(Tp(a + g + 0.5Q));
	}
      return sum;
    }

  template<typename Tp>
    void
    lanczos()
    {
      std::cout.precision(std::numeric_limits<Tp>::digits10);
      std::cout << std::showpoint << std::scientific;
      auto width = 8 + std::cout.precision();

      // From Pugh..
      int n_old = 0;
      int n = -2 - 0.3 * std::log(std::numeric_limits<Tp>::epsilon());
      std::cout << "n_Pugh = " << n << '\n';
      if (n > 1000)
	{
	  std::cerr << "\nlanczos: Calculation of n_Pugh failed.\n";
	  return;
	}

      auto g = n - Tp{0.5Q};
      std::cout << "g = " << g << '\n';
      while (n != n_old)
	{
	  std::cout << '\n';
	  std::vector<Tp> a;
	  for (int k = 1; k <= n; ++k)
	    {
	      for (int j = 1; j <= k; ++j)
        	std::cout << "  C(" << std::setw(2) << k
			    << ", " << std::setw(2) << j
			    << ") = " << std::setw(4) << cheby(k, j);
              std::cout << '\n';
	    }

	  std::cout << '\n';
	  auto prev = std::numeric_limits<Tp>::max();
	  for (int k = 0; k <= n + 5; ++k)
	    {
	      auto curr = p(k, g);
	      if (std::abs(curr) > std::abs(prev))
		{
		  n_old = n;
		  n = k;
		  g = n - Tp{0.5Q};
		  break;
		}
	      prev = curr;
	      std::cout << "  p(" << k << ", " << g << ") = " << curr << '\n';
	    }
	  std::cout << "n = " << n << '\n';
	}

      auto log_gamma1p_lanczos =
	[=](Tp z)
	-> Tp
	{
	  constexpr auto s_ln_2 = emsr::ln2_v<Tp>;
	  constexpr auto s_ln_pi = emsr::lnpi_v<Tp>;
	  constexpr auto s_log_sqrt_2pi = (s_ln_2 + s_ln_pi) / Tp{2};
	  auto fact = Tp{1};
	  auto sum = Tp{0.5Q} * p(0, g);
	  for (int k = 1; k < n; ++k)
	    {
	      fact *= (z - k + 1) / (z + k);
	      sum += fact * p(k, g);
	    }
	  return s_log_sqrt_2pi + std::log(sum)
	       + (z + Tp{0.5Q}) * std::log(z + g + Tp{0.5Q})
	       - (z + g + Tp{0.5Q});
	};

      std::cout << '\n'
		<< ' ' << std::setw(width) << "z"
		<< ' ' << std::setw(width) << "lanczos"
		<< ' ' << std::setw(width) << "lgamma"
		<< ' ' << std::setw(width) << "delta"
		<< '\n';
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = Tp{0.01Q} * i;
	  std::cout << ' ' << std::setw(width) << z
		    << ' ' << std::setw(width) << log_gamma1p_lanczos(z - Tp{1})
		    << ' ' << std::setw(width) << std::lgamma(z)
		    << ' ' << std::setw(width) << log_gamma1p_lanczos(z - Tp{1}) - std::lgamma(z)
		    << '\n';
	}
    }


  /**
   * A struct for Lanczos algorithm Chebyshev arrays of coefficients.
   */
  template<typename Tp>
    class _GammaLanczos
    {
    };

  template<>
    class _GammaLanczos<float>
    {
    public:
    private:
      static constexpr float s_g = 6.5F;
      static constexpr std::array<float, 7>
      s_cheby
      {
	 3.307139e+02F,
	-2.255998e+02F,
	 6.989520e+01F,
	-9.058929e+00F,
	 4.110107e-01F,
	-4.150391e-03F,
	-3.417969e-03F,
      };
    };

  template<>
    class _GammaLanczos<double>
    {
    public:
    private:
      static constexpr double s_g = 9.5;
      static constexpr std::array<double, 10>
      s_cheby
      {
	 5.557569219204146e+03,
	-4.248114953727554e+03,
	 1.881719608233706e+03,
	-4.705537221412237e+02,
	 6.325224688788239e+01,
	-4.206901076213398e+00,
	 1.202512485324405e-01,
	-1.141081476816908e-03,
	 2.055079676210880e-06,
	 1.280568540096283e-09,
      };
    };

  template<>
    class _GammaLanczos<long double>
    {
    public:
    private:
      static constexpr long double s_g = 10.5L;
      static constexpr std::array<long double, 11>
      s_cheby
      {
	 1.440399692024250728e+04L,
	-1.128006201837065341e+04L,
	 5.384108670160999829e+03L,
	-1.536234184127325861e+03L,
	 2.528551924697309561e+02L,
	-2.265389090278717887e+01L,
	 1.006663776178612579e+00L,
	-1.900805731354182626e-02L,
	 1.150508317664389324e-04L,
	-1.208915136885480024e-07L,
	-1.518856151960790157e-10L,
      };
    };

#if !defined(STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  template<>
    class _GammaLanczos<__float128>
    {
    public:
    private:
      static constexpr __float128 s_g = 13.5Q;
      static constexpr std::array<__float128, 14>
      s_cheby
      {
	 2.564476893267270739326759539239521e+05Q,
	-2.115503710351143292058877626137631e+05Q,
	 1.184019145386031178030546551084647e+05Q,
	-4.454725318050137764345651484051269e+04Q,
	 1.108444304734911694218772523380670e+04Q,
	-1.779039702565864879298077142685014e+03Q,
	 1.775917720477127714789670669062174e+02Q,
	-1.046153541164985987873520383376722e+01Q,
	 3.366596234159295518095091028749941e-01Q,
	-5.260241043537793269590387229607226e-03Q,
	 3.290526994052737279641248387319708e-05Q,
	-5.791840871423216610203896349432783e-08Q,
	 1.318192361117106013253744235732287e-11Q,
	-3.352799410216973507737605805778536e-17Q,
      };
    };
#endif

int
main()
{
  std::cout << "\n\nLanczos Algorithm\n\n";
  std::cout << "\nlanczos<float>\n";
  lanczos<float>();
  std::cout << "\nlanczos<double>\n";
  lanczos<double>();
  std::cout << "\nlanczos<long double>\n";
  lanczos<long double>();
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  //FIXME!!! std::cout << "\nlanczos<__float128>\n";
  //FIXME!!! lanczos<__float128>();
#endif
}
