#include <e_float/e_float.h>
#include <examples/examples.h>
#include <functions/functions.h>
#include <utility/util_function_derivative.h>

namespace examples
{
  namespace nr_004
  {
    class BesselJvDerivativeWithRespectTo_v : public Util::FunctionDerivative<e_float>
    {
    private:

      const e_float my_x;

    public:

      BesselJvDerivativeWithRespectTo_v(const e_float& v, const e_float& x) : Util::FunctionDerivative<e_float>(v),
                                                                              my_x(x) { }

      virtual ~BesselJvDerivativeWithRespectTo_v() { }

    private:

      virtual e_float my_function(const e_float& v) const
      {
        return ef::cyl_bessel_j(v, my_x);
      }
    };
  }
}

e_float examples::nr_004::bessel_jv_derivative_wrt_v(const e_float& v, const e_float& x)
{
  const BesselJvDerivativeWithRespectTo_v the_d(v, x);

  return the_d.operation();
}

e_float examples::nr_004::bessel_jv_derivative_wrt_v_test(void)
{
  static const e_float v = 123 + ef::catalan();
  static const e_float x = 151 + ef::euler_gamma();

  // -0.04527629418078693317354422625310156043076923479820609148366251450803148252456061423366518929022315670077714828743749175114999591662663522613406803013423433916263849156834813110673993523679418636346879225816681038135723802070664699170030380250851735369952242936003424802091315905208453804904021597884135812760485879360081229809807866113955718822069367561971830479680406462810714616222515895003673776027
  return bessel_jv_derivative_wrt_v(v, x);
}
