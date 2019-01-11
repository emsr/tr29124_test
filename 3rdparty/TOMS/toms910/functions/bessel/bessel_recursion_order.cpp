
#include <e_float/e_float.h>
#include <functions/bessel/bessel_recursion_order.h>
#include <functions/constants/constants.h>
#include <functions/elementary/elementary.h>
#include <utility/util_find_root_bisect.h>

namespace BesselRecursionOrder
{
  class FindRootStartJx : public Util::FindRootBisect<double>
  {
  private:

    static const double& my_lo (void) { static const double val = static_cast<double>(1.0);   return val; }
    static const double& my_hi (void) { static const double val = static_cast<double>(1.0e7); return val; }
    static const double& my_tol(void) { static const double val = static_cast<double>(0.1);   return val; }

  protected:

    static const double& my_digits(void)
    {
      static const INT32  nd = static_cast<INT32>(ef::tol());
      static const double xd = static_cast<double>(nd);
      return xd;
    }

    FindRootStartJx() : Util::FindRootBisect<double>(my_lo(), my_hi(), my_tol()) { }

  public:
  
    virtual ~FindRootStartJx() { }
  };

  class FindRootStartJ0 : public FindRootStartJx
  {
  protected:

    const double x;

  public:

    FindRootStartJ0(const double& X) : x(X) { }

    virtual ~FindRootStartJ0() { }

  private:

    virtual double ez_term(void) const { return static_cast<double>(1.3591409142295226177) * x; }

    virtual double my_function(const double& v) const
    {
      // This equation corresponds to finding the root value of v for |Jv(z)| = 10^(-digits)
      // using the large order asymptotic behavior of Jv.
      // 
      //                          1             e * z
      // Jv asymptotic --> ----------------  * [-----]^v = 10^-digits
      //                   sqrt(2 * pi * v)     2 * v

      const double log_ez_term_over_v = ::log10(ez_term()) - ::log10(v);

      return     my_digits()
               - (static_cast<double>(0.5) * ::log10(static_cast<double>(6.2831853071795864769) * v))
               + (v * log_ez_term_over_v);
    }
  };

  class FindRootStartJn : public FindRootStartJx
  {
  protected:

    const INT32  n;
    const double x;

  public:

    FindRootStartJn(const INT32 N,
                    const double& X) : n(N < 20 ? 20 : N),
                                       x(X) { }

    virtual ~FindRootStartJn() { }

  private:

    virtual double ez_term(void) const { return static_cast<double>(1.3591409142295226177) * x; }

    virtual double my_function(const double& v) const
    {
      // This equation corresponds to finding the root value of v for |Jv(z)| = 10^(-p / 2) * |Jn(z)|
      // using the large order asymptotic behavior of Jv.
      //
      //                          1             e * z                         1             e * z
      // Jm asymptotic --> ----------------  * [-----]^v = 10^(-p/2) * ----------------  * [-----]^n
      //                   sqrt(2 * pi * v)     2 * v                  sqrt(2 * pi * n)     2 * n
      //

      const double log_ez_term_over_v = ::log10(ez_term()) - ::log10(v);
      const double log_ez_term_over_n = ::log10(ez_term()) - ::log10(static_cast<double>(n));

      return   (static_cast<double>(0.5) * my_digits())
             - (static_cast<double>(0.5) * ::log10(v))
             + (v * log_ez_term_over_v)
             + (static_cast<double>(0.5) * ::log10(static_cast<double>(n)))
             - (static_cast<double>(n) * log_ez_term_over_n);
    }
  };

  class FindRootStartI0 : public FindRootStartJ0
  {
  public:

    FindRootStartI0(const double& X) : FindRootStartJ0(X) { }

    virtual ~FindRootStartI0() { }

  private:

    virtual double ez_term(void) const { return static_cast<double>(2.7182818284590452354) * x; }
  };

  class FindRootStartIn : public FindRootStartJn
  {
  public:

    FindRootStartIn(const INT32 N,
                    const double& X) : FindRootStartJn(N, X) { }

    virtual ~FindRootStartIn() { }

  private:

    virtual double ez_term(void) const { return static_cast<double>(2.7182818284590452354) * x; }
  };

  class FindRootStartKv : public FindRootStartJx
  {
  protected:

    const double v;
    const double z;

  public:

    FindRootStartKv(const double& V,
                    const double& Z) : v(V), z(Z) { }

    virtual ~FindRootStartKv() { }

  private:

    virtual double my_trig_v(void) const { return ::cos(static_cast<double>(3.1415926535897932385) * v); }

    virtual double my_function(const double& n) const
    {
      // This equation corresponds to finding the root value of n for |kn(z)| = 10^(-p)
      // using the large order asymptotic behavior of kn shown in Temme, Equation 3.9.
      // 
      // kn asymptotic --> pi^(-1/2) cos(v pi) 2^(1/4) n^(-1/2) z^(-v - 1/4) exp[z - 2 sqrt(2nz)]
      //

      return   (my_digits() * static_cast<double>(2.3025850929940456840))
             - (static_cast<double>(0.5) * ::log(static_cast<double>(3.1415926535897932385) * n))
             - ((v + static_cast<double>(0.25)) * ::log(z))
             + (static_cast<double>(0.25) * static_cast<double>(0.69314718055994530942))
             + ::log(my_trig_v())
             + z
             - (static_cast<double>(2.0) * ::sqrt(static_cast<double>(2.0) * (n * z)));
    }
  };

  class FindRootStartKv_v_near_half : public FindRootStartKv
  {
  public:

    FindRootStartKv_v_near_half(const double& V,
                                const double& Z) : FindRootStartKv(V, Z) { }

    virtual ~FindRootStartKv_v_near_half() { }

  private:

    // This equation corresponds to the equation which is defined for RecursionKv.
    // However, an approximation of cos(v pi) --> sin(v) for v very close to 1/2.

    virtual double my_trig_v(void) const { return ::sin(v); }
  };
}

INT32 BesselRecursionOrder::RecursionStartOrderJ0(const double x)
{
  const FindRootStartJ0 find_root(x);
  const double start_order = find_root.operation();
  return static_cast<INT32>(static_cast<INT64>(start_order)) + static_cast<INT32>(20);
}

INT32 BesselRecursionOrder::RecursionStartOrderJn(const INT32 n, const double x)
{
  const FindRootStartJn find_root(n, x);
  const double start_order = find_root.operation();
  return static_cast<INT32>(static_cast<INT64>(start_order)) + static_cast<INT32>(20);
}

INT32 BesselRecursionOrder::RecursionStartOrderI0(const double x)
{
  const FindRootStartI0 find_root(x);
  const double start_order = find_root.operation();
  return static_cast<INT32>(static_cast<INT64>(start_order)) + static_cast<INT32>(20);
}

INT32 BesselRecursionOrder::RecursionStartOrderIn(const INT32 n, const double x)
{
  const FindRootStartIn find_root(n, x);
  const double start_order = find_root.operation();
  return static_cast<INT32>(static_cast<INT64>(start_order)) + static_cast<INT32>(20);
}

INT32 BesselRecursionOrder::RecursionStartOrderKv(const double v, const double x)
{
  const FindRootStartKv find_root(v, x);
  const double start_order = find_root.operation();
  return static_cast<INT32>(static_cast<INT64>(start_order)) + static_cast<INT32>(20);
}

INT32 BesselRecursionOrder::RecursionStartOrderKv_v_near_half(const double v, const double x)
{
  const FindRootStartKv_v_near_half find_root(v, x);
  const double start_order = find_root.operation();
  return static_cast<INT32>(static_cast<INT64>(start_order)) + static_cast<INT32>(20);
}
