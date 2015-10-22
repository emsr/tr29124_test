
// ../bin/bin/g++ -std=c++0x -o test_complex_gamma test_complex_gamma.cpp

//  This works.  We could replace the lanczos in tr1 with this one.
// You would define clgamma[f,,l], ctgamma[f,,l] also.

//  Strange, the __numeric_constants work with complex template args!
//  OK...
//      static _Tp __pi() throw()
//      { return static_cast<_Tp>(3.1415926535897932384626433832795029L); }

//  Check this!!!!
//Notes
//POSIX specification additionally requires that each execution of lgamma
//stores the sign of the gamma function of arg in the external variable signgam.
//  Should i do that shit?

// http://en.wikipedia.org/wiki/Particular_values_of_the_gamma_function

#include <complex>
#include <iostream>
#include <stdexcept>
#include <tr1/cmath>

using namespace std::tr1::__detail;

template<typename _Tp>
  struct __value_type
  {
    typedef _Tp value_type;
  };

template<typename _Tp>
  struct __value_type<std::complex<_Tp>>
  {
    typedef typename std::complex<_Tp>::value_type value_type;
  };

// keep the sign for complex arguments...
template<typename _Tp>
  _Tp __sin_fact_helper(_Tp __x)
  { return std::abs(__x); };

template<typename _Tp>
  std::complex<_Tp> __sin_fact_helper(std::complex<_Tp> __x)
  { return __x; };

//  Offer this for real valued lgamma?!?!?:
template<typename _Tp>
  _Tp
  __sign_gamma(_Tp __x)
  {
    if (__x < _Tp(0) && int(__x) % 2 == 0)
      return -_Tp(1);
    else
      return _Tp(1);
  }

    /**
     *   @brief Return \f$log(\Gamma(x))\f$ by the Lanczos method.
     *          This method dominates all others on the positive axis I think.
     *
     *   @param __x The argument of the log of the gamma function.
     *   @return  The logarithm of the gamma function.
     */
    template<typename _Tp>
    _Tp
    __log_cccgamma_lanczos(_Tp __x)
    {
      typedef typename __value_type<_Tp>::value_type _Val;
      const _Tp __xm1 = __x - _Tp(1);

      static const _Val __lanczos_cheb_7[9] = {
       _Val( 0.99999999999980993227684700473478L),
       _Val( 676.520368121885098567009190444019L),
       _Val(-1259.13921672240287047156078755283L),
       _Val( 771.3234287776530788486528258894L),
       _Val(-176.61502916214059906584551354L),
       _Val( 12.507343278686904814458936853L),
       _Val(-0.13857109526572011689554707L),
       _Val( 9.984369578019570859563e-6L),
       _Val( 1.50563273514931155834e-7L)
      };

      static const _Val __LOGROOT2PI
          = _Val(0.9189385332046727417803297364056176L);

      _Tp __sum = __lanczos_cheb_7[0];
      for(unsigned int __k = 1; __k < 9; ++__k)
        __sum += __lanczos_cheb_7[__k] / (__xm1 + _Tp(__k));

      const _Tp __term1 = (__xm1 + _Tp(0.5L))
                        * std::log((__xm1 + _Tp(7.5L))
                       / __numeric_constants<_Tp>::__euler());
      const _Tp __term2 = __LOGROOT2PI + std::log(__sum);
      const _Tp __result = __term1 + (__term2 - _Tp(7L));

      return __result;
    }


    /**
     *   @brief Return \f$ log(|\Gamma(x)|) \f$.
     *          This will return values even for \f$ x < 0 \f$.
     *          To recover the sign of \f$ \Gamma(x) \f$ for
     *          any argument use @a __log_gamma_sign.
     *
     *   @param __x The argument of the log of the gamma function.
     *   @return  The logarithm of the gamma function.
     */
    template<typename _Tp>
    _Tp
    __log_cccgamma(_Tp __x)
    {
      typedef typename __value_type<_Tp>::value_type _Val;
      if (std::real(__x) > _Val(0.5L))
        return __log_cccgamma_lanczos(__x);
      else
        {
          const _Tp __sin_fact
                 = __sin_fact_helper(std::sin(__numeric_constants<_Tp>::__pi() * __x));
          if (__sin_fact == _Tp(0))
            //std::__throw_domain_error(__N("Argument is nonpositive integer "
            //                              "in __log_gamma"));
            throw std::domain_error("Argument is nonpositive integer in __log_gamma");
          return __numeric_constants<_Tp>::__lnpi()
                     - std::log(__sin_fact)
                     - __log_cccgamma_lanczos(_Tp(1) - __x);
        }
    }

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::complex<double> z(-1.5, 0.5);
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = 1.0/3.0;
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = std::complex<double>(1.0,1.0);
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = std::complex<double>(1.0,-1.0);
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = std::complex<double>(0.5,0.5);
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = std::complex<double>(0.5,-0.5);
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = std::complex<double>(5.0,3.0);
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = std::complex<double>(5.0,-3.0);
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;

  //  Test reals...
  double x = 1.0/3.0;
  std::cout << "log_gamma(" << x << ") = " << std::exp(__log_cccgamma(x)) << std::endl;
  x = 1.0/4.0;
  std::cout << "log_gamma(" << x << ") = " << std::exp(__log_cccgamma(x)) << std::endl;
  x = 1.0/5.0;
  std::cout << "log_gamma(" << x << ") = " << std::exp(__log_cccgamma(x)) << std::endl;
  x = 1.0/6.0;
  std::cout << "log_gamma(" << x << ") = " << std::exp(__log_cccgamma(x)) << std::endl;
  x = 1.0/7.0;
  std::cout << "log_gamma(" << x << ") = " << std::exp(__log_cccgamma(x)) << std::endl;

  //  Test sign with complex...  got it with the helper!  No need for a sign_gamma() or a signgam global variable.
  z = -0.5;
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = -1.5;
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = -2.5;
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = -3.5;
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;
  z = -4.5;
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;

  //  Test sign_gamma()...
  std::cout << "__sign_gamma(" << -0.5 << ") = " << __sign_gamma(-0.5) << std::endl;
  std::cout << "__sign_gamma(" << -1.5 << ") = " << __sign_gamma(-1.5) << std::endl;
  std::cout << "__sign_gamma(" << -2.5 << ") = " << __sign_gamma(-2.5) << std::endl;
  std::cout << "__sign_gamma(" << -3.5 << ") = " << __sign_gamma(-3.5) << std::endl;
  std::cout << "__sign_gamma(" << -4.5 << ") = " << __sign_gamma(-4.5) << std::endl;

  //  Test reflection...
  z = std::complex<double>(-8.5, 0.5);
  std::cout << "log_gamma(" << z << ") = " << std::exp(__log_cccgamma(z)) << std::endl;

  //  WTF...
  std::cout << "pi      = " << __numeric_constants<std::complex<double>>::__pi() << std::endl;
  std::cout << "gamma_E = " << __numeric_constants<std::complex<double>>::__euler() << std::endl;
  std::cout << "ln(pi)  = " << __numeric_constants<std::complex<double>>::__lnpi() << std::endl;
}
