#include <functions/complex/e_float_complex.h>
#include <examples/examples.h>
#include <functions/constants/constants.h>
//#include <interop/cas/mathematica/ef_cas_object_mathematica.h>
#include <interop/cas/dummy/ef_cas_object_dummy.h>
#include <utility/util_lexical_cast.h>

namespace examples
{
  namespace nr_007
  {
    template<typename T> static inline T poly_gamma_mathematica_template(const T& v, const T& z)
    {
      const ef_cas::ComputerAlgebraSystemObject& the_cas = ef_cas::cas_dummy();
//      const ef_cas::ComputerAlgebraSystemObject& the_cas = ef_cas::cas_mathematica();

      std::string str_v;
      std::string str_z;

      the_cas.create_mp_string(str_v, v);
      the_cas.create_mp_string(str_z, z);

      // Create a string such as: Table[N[PolyLog[((Re_v`115) + ((Im_v`115) I)), ((Re_z`115) + ((Im_z`115) I))], 130], {1}]
      const std::string str_cmd = "Table[N[PolyGamma[" + str_v + ", " + str_z + "], " + Util::lexical_cast(the_cas.my_digits()) + "], {1}]";

      std::vector<T> values;
      const bool b_get_value = the_cas.get_values(str_cmd, values);

      return (b_get_value ? values.front() : std::numeric_limits<e_float>::quiet_NaN());
    }
  }
}

e_float examples::nr_007::poly_gamma_mathematica(const e_float& v, const e_float& x)
{
  return poly_gamma_mathematica_template(v, x);
}

ef_complex examples::nr_007::poly_gamma_mathematica(const ef_complex& v, const ef_complex& z)
{
  return poly_gamma_mathematica_template(v, z);
}

ef_complex examples::nr_007::poly_gamma_mathematica_test(void)
{
  static const ef_complex v(ef::euler_gamma(), ef::third());
  static const ef_complex z(e_float(123) / 7, ef::catalan());

  return poly_gamma_mathematica(v, z);
}
