#ifndef SF_HYPERG_H
#define SF_HYPERG_H 1

#include <emsr/fp_type_util.h>
#include <emsr/detail/sf_hyperg.tcc>

namespace emsr
{

  // Confluent hypergeometric functions

  /**
   * Return the confluent hypergeometric function @f$ {}_1F_1(a;c;x) @f$
   * of real numerator parameter @c a, denominator parameter @c c,
   * and argument @c x.
   *
   * The confluent hypergeometric function is defined by
   * @f[
   *    {}_1F_1(a;c;x) = \sum_{n=0}^{\infty} \frac{(a)_n x^n}{(c)_n n!}
   * @f]
   * where the Pochhammer symbol is @f$ (x)_k = (x)(x+1)...(x+k-1) @f$,
   * @f$ (x)_0 = 1 @f$
   *
   * @param a The numerator parameter
   * @param c The denominator parameter
   * @param x The argument
   */
  template<typename Tpa, typename Tpc, typename Tp>
    inline emsr::fp_promote_t<Tpa, Tpc, Tp>
    conf_hyperg(Tpa a, Tpc c, Tp x)
    {
      using type = emsr::fp_promote_t<Tpa, Tpc, Tp>;
      return emsr::detail::conf_hyperg<type>(a, c, x);
    }

  // Confluent hypergeometric functions

  /**
   * Return the Tricomi confluent hypergeometric function @f$ U(a,c,x) @f$
   * of real numerator parameter @c a, denominator parameter @c c,
   * and argument @c x.
   *
   * The Tricomi confluent hypergeometric function is defined by
   * @f[
   *    U(a,c,x) = \frac{\Gamma(1-c)}{\Gamma(a-c+1)} {}_1F_1(a;c;x)
   *       + \frac{\Gamma(c-1)}{\Gamma(a)} x^{1-c} {}_1F_1(a-c+1;2-c;x)
   * @f]
   * where @f$ {}_1F_1(a;c;x) @f$ if the confluent hypergeometric function.
   *
   * @see conf_hyperg.
   *
   * @param a The numerator parameter
   * @param c The denominator parameter
   * @param x The argument
   */
  template<typename Tpa, typename Tpc, typename Tp>
    inline emsr::fp_promote_t<Tpa, Tpc, Tp>
    tricomi_u(Tpa a, Tpc c, Tp x)
    {
      using type = emsr::fp_promote_t<Tpa, Tpc, Tp>;
      return emsr::detail::tricomi_u<type>(a, c, x);
    }

  // Hypergeometric functions

  /**
   * Return the hypergeometric function @f$ {}_2F_1(a,b;c;x) @f$
   * of real numerator parameters @c a and @c b,
   * denominator parameter @c c, and argument @c x.
   *
   * The hypergeometric function is defined by
   * @f[
   *    {}_2F_1(a,b;c;x) = \sum_{n=0}^{\infty} \frac{(a)_n (b)_n x^n}{(c)_n n!}
   * @f]
   * where the Pochhammer symbol is @f$ (x)_k = (x)(x+1)...(x+k-1) @f$,
   * @f$ (x)_0 = 1 @f$
   *
   * @param a The first numerator parameter
   * @param b The second numerator parameter
   * @param c The denominator parameter
   * @param x The argument
   */
  template<typename Tpa, typename Tpb, typename Tpc, typename Tp>
    inline typename emsr::fp_promote_t<Tpa, Tpb, Tpc, Tp>
    hyperg(Tpa a, Tpb b, Tpc c, Tp x)
    {
      using type = emsr::fp_promote_t<Tpa, Tpb, Tpc, Tp>;
      return emsr::detail::hyperg<type>(a, b, c, x);
    }

  // Confluent hypergeometric limit functions

  /**
   * Return the confluent hypergeometric limit function @f$ {}_0F_1(;c;x) @f$
   * of real numerator parameter @c c and argument @c x.
   *
   * The confluent hypergeometric limit function is defined by
   * @f[
   *    {}_0F_1(;c;x) = \sum_{n=0}^{\infty} \frac{x^n}{(c)_n n!}
   * @f]
   * where the Pochhammer symbol is @f$ (x)_k = (x)(x+1)...(x+k-1) @f$,
   * @f$ (x)_0 = 1 @f$
   *
   * @param c The denominator parameter
   * @param x The argument
   */
  template<typename Tpc, typename Tp>
    inline typename emsr::fp_promote_t<Tpc, Tp>
    conf_hyperg_lim(Tpc c, Tp x)
    {
      typedef typename emsr::fp_promote_t<Tpc, Tp> type;
      return emsr::detail::conf_hyperg_lim<type>(c, x);
    }

} // namespace emsr

#endif // SF_HYPERG_H
