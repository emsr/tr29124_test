#include <e_float/e_float.h>
#include <examples/examples.h>
#include <functions/functions.h>
#include <utility/util_lexical_cast.h>
#include <utility/util_trapezoid.h>

namespace examples
{
  namespace nr_005
  {
    class RecursiveTrapezoidJ0 : public Util::RecursiveTrapezoidRule<e_float>
    {
    private:

      static const e_float& my_tol(void)
      {
        static const e_float val("1E-" + Util::lexical_cast(std::numeric_limits<e_float>::digits10 / 2));
        return val;
      }

    private:

      const e_float my_z;

    public:

      RecursiveTrapezoidJ0(const e_float& z) : Util::RecursiveTrapezoidRule<e_float>(ef::zero(), ef::pi(), my_tol()),
                                               my_z(z) { }

      virtual ~RecursiveTrapezoidJ0() { }

    private:

      virtual e_float my_function(const e_float& x) const
      {
        return ef::cos(my_z * ef::sin(x));
      }
    };
  }
}

e_float examples::nr_005::recursive_trapezoid_j0(const e_float& x)
{
  const RecursiveTrapezoidJ0 rtj0(x);

  return rtj0.operation() / ef::pi();
}

e_float examples::nr_005::recursive_trapezoid_j0_test(void)
{
  static const e_float x = 12 + ef::euler_gamma();

  return recursive_trapezoid_j0(x);
}
