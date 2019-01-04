
#ifndef _BESSEL_ASYMP_2008_04_09_H_
  #define _BESSEL_ASYMP_2008_04_09_H_

  #include <vector>

  #include <e_float/e_float.h>
  #include <functions/constants/constants.h>
  #include <functions/elementary/elementary.h>
  #include <functions/bessel/airy.h>
  #include <utility/util_digit_scale.h>

  namespace BesselAsymp
  {
    struct Method
    {
    private:

      // Non-copyable
      const Method& operator=(const Method&);
      Method(const Method&);

    protected:

      Method() { }

    public:

      virtual ~Method() { }

      virtual INT32 lower(const double& xd) const = 0; // NOCOVER_LINE
      virtual INT32 upper(const double& xd) const = 0; // NOCOVER_LINE

      virtual e_float calc(const e_float& v, const e_float& x) const = 0; // NOCOVER_LINE

      bool within(const INT32 n, const double x) const
      {
        const INT32 n_lo = lower(x);
        const INT32 n_hi = upper(x);

        return (n >= n_lo) && (n <= n_hi);
      }

      static double xmin(void)
      {
        return static_cast<double>(8000.0) * Util::DigitScale();
      }
    };

    struct UniformAsyBase : public Method
    {
      UniformAsyBase() { }
      virtual ~UniformAsyBase() { }

      e_float UniformAsymptotic(const e_float& v, const e_float& x) const;

      virtual e_float calc(const e_float& v, const e_float& x) const { return UniformAsymptotic(v, x); }

      virtual bool z_is_gt_one(void) const { return false; }

      virtual e_float my_airy      (const e_float& x) const { return ef::airy_a(x); }
      virtual e_float my_airy_prime(const e_float& x) const { return ef::airy_a_prime(x); }

      virtual bool is_bess_y(void) const { return false; }
    };

    struct UniformAsy_z_gt_one : public UniformAsyBase
    {
      UniformAsy_z_gt_one() { }
      virtual ~UniformAsy_z_gt_one() { }

      virtual bool z_is_gt_one(void) const { return true; }

      virtual INT32 lower(const double& xd) const;
      virtual INT32 upper(const double& xd) const;
    };

    struct UniformAsy_z_lt_one : public UniformAsyBase
    {
      UniformAsy_z_lt_one() { }
      virtual ~UniformAsy_z_lt_one() { }

      virtual INT32 lower(const double& xd) const;
      virtual INT32 upper(const double& xd) const;
    };

    struct DebyeAsyBase : public Method
    {
      DebyeAsyBase() { }
      virtual ~DebyeAsyBase() { }

      e_float DebyeAsymptotic(const e_float& v, const e_float& x) const;

      virtual e_float calc(const e_float& v, const e_float& x) const { return DebyeAsymptotic(v, x); }

      virtual bool is_bess_y(void) const { return false; }

      virtual const e_float& pi_term(void) const { return ef::two_pi(); }
    };

    struct DebyeAsy_BessJ : public DebyeAsyBase
    {
      DebyeAsy_BessJ() { }
      virtual ~DebyeAsy_BessJ() { }
    };

    struct DebyeAsy_BessY : public DebyeAsyBase
    {
      DebyeAsy_BessY() { }
      virtual ~DebyeAsy_BessY() { }

      virtual bool is_bess_y(void) const { return true; }

      virtual const e_float& pi_term(void) const { return ef::pi_half(); }
    };
  }

  namespace BesselJ_AsympMethod
  {
    struct UniformAsy_z_gt_one : public BesselAsymp::UniformAsy_z_gt_one
    {
      UniformAsy_z_gt_one() { }
      virtual ~UniformAsy_z_gt_one() { }
    };

    struct UniformAsy_z_lt_one : public BesselAsymp::UniformAsy_z_lt_one
    {
      UniformAsy_z_lt_one() { }
      virtual ~UniformAsy_z_lt_one() { }
    };
  }

  namespace BesselY_AsympMethod
  {
    struct UniformAsy_z_gt_one : public BesselAsymp::UniformAsy_z_gt_one
    {
      UniformAsy_z_gt_one() { }
      virtual ~UniformAsy_z_gt_one() { }

      virtual e_float my_airy      (const e_float& x) const { return ef::airy_b(x); }
      virtual e_float my_airy_prime(const e_float& x) const { return ef::airy_b_prime(x); }

      virtual bool is_bess_y(void) const { return true; }
    };

    struct UniformAsy_z_lt_one : public BesselAsymp::UniformAsy_z_lt_one
    {
      UniformAsy_z_lt_one() { }
      virtual ~UniformAsy_z_lt_one() { }

      virtual e_float my_airy      (const e_float& x) const { return ef::airy_b(x); }
      virtual e_float my_airy_prime(const e_float& x) const { return ef::airy_b_prime(x); }

      virtual bool is_bess_y(void) const { return true; }
    };
  }

  namespace BesselAsymp
  {
    const std::vector<e_float>& LambdaS(void);
    const std::vector<e_float>& MuS    (void);

    e_float DebyeU(const INT32 n,
                   const e_float& t,
                   std::deque<e_float>& vt,
                   const bool has_phase = false);
  }

#endif // _BESSEL_ASYMP_2008_04_09_H_

