/**
 *
 */

//  This works.  We could replace the lanczos in tr1 with this one.
// You would define clgamma[f,,l], ctgamma[f,,l] also.

//  Check this!!!!
//Notes
//POSIX specification additionally requires that each execution of lgamma
//stores the sign of the gamma function of arg in the external variable signgam.
//  Should i do that shit?

// http://en.wikipedia.org/wiki/Particular_values_of_the_gamma_function

#include <complex>
#include <iostream>
#include <stdexcept>
#include <cmath>
#include <limits>

#include <emsr/math_constants.h>

template<typename Tp>
  struct value_type
  {
    using type = Tp;
  };

template<typename Tp>
  struct value_type<std::complex<Tp>>
  {
    using type = typename std::complex<Tp>::value_type;
  };

// keep the sign for complex arguments...
template<typename Tp>
  Tp sin_fact_helper(Tp x)
  { return std::abs(x); };

template<typename Tp>
  std::complex<Tp> sin_fact_helper(std::complex<Tp> x)
  { return x; };

//  Offer this for real valued lgamma?!?!?:
template<typename Tp>
  Tp
  sign_gamma(Tp x)
  {
    if (x < Tp(0) && int(x) % 2 == 0)
      return -Tp(1);
    else
      return Tp(1);
  }

    /**
     *   @brief Return \f$log(\Gamma(x))\f$ by the Lanczos method.
     *          This method dominates all others on the positive axis I think.
     *
     *   @param x The argument of the log of the gamma function.
     *   @return  The logarithm of the gamma function.
     */
    template<typename Tp>
    Tp
    log_cccgamma_lanczos(Tp x)
    {
      typedef typename value_type<Tp>::type Val;
      const Tp xm1 = x - Tp(1);

      static const Val lanczos_cheb_7[9] = {
       Val( 0.99999999999980993227684700473478L),
       Val( 676.520368121885098567009190444019L),
       Val(-1259.13921672240287047156078755283L),
       Val( 771.3234287776530788486528258894L),
       Val(-176.61502916214059906584551354L),
       Val( 12.507343278686904814458936853L),
       Val(-0.13857109526572011689554707L),
       Val( 9.984369578019570859563e-6L),
       Val( 1.50563273514931155834e-7L)
      };

      static const Val LOGROOT2PI
          = Val(0.9189385332046727417803297364056176L);

      Tp sum = lanczos_cheb_7[0];
      for(unsigned int k = 1; k < 9; ++k)
        sum += lanczos_cheb_7[k] / (xm1 + Tp(k));

      const Tp term1 = (xm1 + Tp(0.5L))
                        * std::log((xm1 + Tp(7.5L))
                       / emsr::egamma_v<Val>);
      const Tp term2 = LOGROOT2PI + std::log(sum);
      const Tp result = term1 + (term2 - Tp(7L));

      return result;
    }


    /**
     *   @brief Return \f$ log(|\Gamma(x)|) \f$.
     *          This will return values even for \f$ x < 0 \f$.
     *          To recover the sign of \f$ \Gamma(x) \f$ for
     *          any argument use @a log_gamma_sign.
     *
     *   @param x The argument of the log of the gamma function.
     *   @return  The logarithm of the gamma function.
     */
    template<typename Tp>
    Tp
    log_cccgamma(Tp x)
    {
      typedef typename value_type<Tp>::type Val;
      if (std::real(x) > Val(0.5L))
        return log_cccgamma_lanczos(x);
      else
        {
          const Tp sin_fact
                 = sin_fact_helper(std::sin(emsr::pi_v<Val> * x));
          if (sin_fact == Tp(0))
            //throw std::domain_error("Argument is nonpositive integer in log_gamma");
            throw std::domain_error("Argument is nonpositive integer in log_gamma");
          return emsr::lnpi_v<Val>
                     - std::log(sin_fact)
                     - log_cccgamma_lanczos(Tp(1) - x);
        }
    }

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::complex<double> z(-1.5, 0.5);
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = 1.0/3.0;
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = std::complex<double>(1.0,1.0);
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = std::complex<double>(1.0,-1.0);
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = std::complex<double>(0.5,0.5);
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = std::complex<double>(0.5,-0.5);
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = std::complex<double>(5.0,3.0);
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = std::complex<double>(5.0,-3.0);
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';

  //  Test reals...
  double x = 1.0/3.0;
  std::cout << "gamma(" << x << ") = " << std::exp(log_cccgamma(x)) << '\n';
  x = 1.0/4.0;
  std::cout << "gamma(" << x << ") = " << std::exp(log_cccgamma(x)) << '\n';
  x = 1.0/5.0;
  std::cout << "gamma(" << x << ") = " << std::exp(log_cccgamma(x)) << '\n';
  x = 1.0/6.0;
  std::cout << "gamma(" << x << ") = " << std::exp(log_cccgamma(x)) << '\n';
  x = 1.0/7.0;
  std::cout << "gamma(" << x << ") = " << std::exp(log_cccgamma(x)) << '\n';

  //  Test sign with complex...  got it with the helper!  No need for a sign_gamma() or a signgam global variable.
  z = -0.5;
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = -1.5;
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = -2.5;
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = -3.5;
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
  z = -4.5;
  std::cout << "gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';

  //  Test sign_gamma()...
  std::cout << "sign_gamma(" << -0.5 << ") = " << sign_gamma(-0.5) << '\n';
  std::cout << "sign_gamma(" << -1.5 << ") = " << sign_gamma(-1.5) << '\n';
  std::cout << "sign_gamma(" << -2.5 << ") = " << sign_gamma(-2.5) << '\n';
  std::cout << "sign_gamma(" << -3.5 << ") = " << sign_gamma(-3.5) << '\n';
  std::cout << "sign_gamma(" << -4.5 << ") = " << sign_gamma(-4.5) << '\n';

  //  Test reflection...
  z = std::complex<double>(-8.5, 0.5);
  std::cout << "log_gamma(" << z << ") = " << std::exp(log_cccgamma(z)) << '\n';
}
