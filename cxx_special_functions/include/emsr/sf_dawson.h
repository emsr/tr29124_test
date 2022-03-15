#ifndef SF_DAWSON_H
#define SF_DAWSON_H 1


#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_dawson.tcc>

namespace emsr
{

  // Dawson integral

  /**
   * Return the Dawson integral, @f$ F(x) @f$, for real argument @c x.
   *
   * The Dawson integral is defined by:
   * @f[
   *    F(x) = e^{-x^2}\int_0^x e^{y^2}dy
   * @f]
   * and it's derivative is:
   * @f[
   *    F'(x) = 1 - 2xF(x)
   * @f]
   *
   * @param x The argument @f$ -inf < x < inf @f$.
   */
  template<typename Tp>
    inline emsr::fp_promote_t<Tp>
    dawson(Tp x)
    {
      using type = emsr::fp_promote_t<Tp>;
      return emsr::detail::dawson<type>(x);
    }

} // namespace emsr

#endif // SF_DAWSON_H
