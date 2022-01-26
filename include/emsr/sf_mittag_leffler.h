#ifndef SF_MITTAG_LEFFLER_H
#define SF_MITTAG_LEFFLER_H 1

#include <complex>

#include <emsr/fp_type_util.h>
#include <emsr/sf_mittag_leffler.tcc>

namespace emsr
{

  /**
   * Compute the Mittag-Leffer function:
   * @f[
   *   E_{\alpha,\beta}(z)
   *     = \sum_{k=0}^\infty \frac{z^k}{\Gamma(\beta + \alpha k)},
   *   \mbox{  } \alpha > 0, \beta \in \mathbb{C}, z \in \mathbb{C}
   * @f]
   *
   * @see COMPUTATION OF THE MITTAG-LEFFLER FUNCTION @f$ E_{\alpha,\beta}(z) @f$
   * AND ITS DERIVATIVE, Rudolf Gorenflo, Joulia Loutchko & Yuri Luchko
   */
  template<typename _Tp, typename _Ta, typename _Tb>
    inline std::complex<emsr::fp_promote_t<_Tp, _Ta, _Tb>>
    mittag_leffler(_Ta alpha, _Tb beta, const std::complex<_Tp>& z)
    {
      using type = emsr::fp_promote_t<_Tp, _Ta, _Tb>;
      return emsr::detail::mittag_leffler<type>(alpha, beta, z);
    }

} // namespace emsr

#endif // SF_MITTAG_LEFFLER_H
