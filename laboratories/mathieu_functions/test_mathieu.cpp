/**
 * Accurate Computation of Mathieu Functions,
 * Malcolm M. Bibby and Andrew F. Peterson
 * Morgan & Claypool, 2014
 */

  /**
   * Return the Mathieu function eigenvalue @f$  @f$ .
   */
  template<typename _Tp>
    _Tp
    mathieu_a(int n, _Tp q);

  /**
   * Return the Mathieu function eigenvalue @f$  @f$ .
   */
  template<typename _Tp>
    _Tp
    mathieu_b(int n, _Tp q);

  /**
   * Return the angular Mathieu function @f$ se_n(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_se(int n, _Tp q, _Tp x);

  /**
   * Return the angular Mathieu function @f$ ce_n(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_fe(int n, _Tp q, _Tp x);

  /**
   * Return the angular Mathieu function @f$ fe_n(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_ge(int n, _Tp q, _Tp x);

  /**
   * Return the angular Mathieu function @f$ ge_n(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_ce(int n, _Tp q, _Tp x);

  /**
   * Return the radial Mathieu function @f$ Ms_n^{(k)}(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_ms(int n, int k, _Tp q, _Tp x);

  /**
   * Return the radial Mathieu function @f$ Mc_n^{(k)}(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_mc(int n, int k, _Tp q, _Tp x);

int
main()
{
  return 0;
}
