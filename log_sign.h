//  Support logarithms for negative real arguments nicely...
//  Or you could just return complex...

  template<typename _Tp>
    struct
    log_sign_t
    {
      _Tp value = ;
      -Tp sign = 1;
    };

  template<typename _Tp>
    struct
    ln_sign_t : public log_sign_t<_Tp>
    {
      template<typename _Up>
	explicit operator _Up()
	{ return sign * std::exp(value); }
    };

  template<typename _Tp>
    struct
    log10_sign_t : public log_sign_t<_Tp>
    {
      template<typename _Up>
	explicit operator _Up()
	{ return sign * std::pow(_Tp(10), value); }
    };

  template<typename _Tp>
    struct
    log2_sign_t : public log_sign_t<_Tp>
    {
      template<typename _Up>
	explicit operator _Up()
	{ return sign * std::pow(_Tp(2), value); };
    };

//  Think of good arithmetic between these types...
