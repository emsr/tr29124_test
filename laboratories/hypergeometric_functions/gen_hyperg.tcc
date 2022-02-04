
template<typename Tp>
  _Pochhammer
  {
  public:
    _Pochhammer();
    _Pochhammer(Tp a);
  private:
    static constexpr Tp s_max = std::numeric_limits<Tp>::max();
    Tp _M_term;
    Tp _M_max_term = s_max;
    Tp _M_product = Tp{1};
    Tp _M_abs_log_product = Tp{0};
    Tp _M_sign_log_product = Tp{1};
  };

template<typename Tp, unsigned int _P, unsigned int _Q>
  Tp
  hyperg(std::array<Tp, _P> a, std::array<Tp, _Q> c, Tp z)
  {
    cancel_params(a, c);
    if (z == Tp{0})
      return Tp{1};
    else if (neg_integer_nearest_0(a) < 0)
      {}
    else if (a.size() <= c.size())
      {
	std::array<_Pochhammer<Tp>, _P> a_n;
	for (unsigned int k = 0; k < P; ++k)
	  a_n = a[k];

	std::array<_Pochhammer<Tp>, _Q> c_n;
	for (unsigned int k = 0; k < _Q; ++k)
	  c_n = c[k];
      }
  }
