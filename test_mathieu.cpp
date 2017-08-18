
  /**
   * Return the Mathieu function eigenvalue @f$  @f$ .
   */
  template<typename _Tp>
    _Tp
    mathieu_a(int __n, _Tp __q);

  /**
   * Return the Mathieu function eigenvalue @f$  @f$ .
   */
  template<typename _Tp>
    _Tp
    mathieu_b(int __n, _Tp __q);

  /**
   * Return the angular Mathieu function @f$ se_n^{(k)}(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_se(int __n, int __k, _Tp __q, _Tp __x);

  /**
   * Return the angular Mathieu function @f$ ce_n^{(k)}(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_ce(int __n, int __k, _Tp __q, _Tp __x);

  /**
   * Return the radial Mathieu function @f$ Ms_n^{(k)}(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_ms(int __n, int __k, _Tp __q, _Tp __x);

  /**
   * Return the radial Mathieu function @f$ Mc_n^{(k)}(q,x) @f$.
   */
  template<typename _Tp>
    _Tp
    mathieu_mc(int __n, int __k, _Tp __q, _Tp __x);
