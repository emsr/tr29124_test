

  template <typename _Tp>
  _Tp
  __legendre_p(const int __l, const _Tp __x)
  {

    if (__l < 0)
      throw std::domain_error("Bad order in legendre_p.");
    if (__x < -_Tp(1) || __x > +_Tp(1))
      throw std::domain_error("Bad argument x in legendre_p.");

    _Tp __p_lm2 = _Tp(1);
    if (__l == _Tp(0))
      return __p_lm2;

    _Tp __p_lm1 = __x;
    if (__l == _Tp(1))
      return __p_lm1;

    if (__x == +_Tp(1))
      return +_Tp(1);

    if (__x == -_Tp(1))
      return (__l % 2 == 1 ? -_Tp(1) : +_Tp(1));

    _Tp __p_l;
    for (int __ll = 2; __ll <= __l; ++__ll)
      {
        __p_l = (_Tp(2 * __ll - 1) * __x * __p_lm1
               - _Tp(__ll - 1) * __p_lm2) / _Tp(__ll);
        __p_lm2 = __p_lm1;
        __p_lm1 = __p_l;
      }

    return __p_l;
  }
