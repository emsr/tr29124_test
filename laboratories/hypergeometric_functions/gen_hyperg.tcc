
template<typename _Tp>
  _Pochhammer
  {
  public:
    _Pochhammer();
    _Pochhammer(_Tp a);
  private:
    static constexpr _Tp s_max = std::numeric_limits<_Tp>::max();
    _Tp _M_term;
    _Tp _M_max_term = s_max;
    _Tp _M_product = _Tp{1};
    _Tp _M_abs_log_product = _Tp{0};
    _Tp _M_sign_log_product = _Tp{1};
  };

template<typename _Tp, unsigned int _P, unsigned int _Q>
  Tp
  hyperg(std::array<_Tp, _P> a, std::array<_Tp, _Q> c, _Tp z)
  {
    cancel_params(a, c);
    if (z == _Tp{0})
      return _Tp{1};
    else if (neg_integer_nearest_0(a) < 0)
      {}
    else if (a.size() <= c.size())
      {
	std::array<_Pochhammer<_Tp>, _P> a_n;
	for (unsigned int k = 0; k < P; ++k)
	  a_n = a[k];

	std::array<_Pochhammer<_Tp>, _Q> c_n;
	for (unsigned int k = 0; k < _Q; ++k)
	  c_n = c[k];
      }
  }
