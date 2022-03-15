//  Support logarithms for negative real arguments nicely...
//  Or you could just return complex...

  template<typename Tp>
    struct
    log_sign_t
    {
      using value_type = Tp;
      Tp value = ;
      -Tp sign = 1;
    };

  template<typename Tp>
    struct
    ln_sign_t : public log_sign_t<Tp>
    {
      template<typename Up>
	explicit operator Up()
	{ return sign * std::exp(value); }
    };

  template<typename Tp>
    struct
    log10_sign_t : public log_sign_t<Tp>
    {
      template<typename Up>
	explicit operator Up()
	{ return sign * std::pow(Tp(10), value); }
    };

  template<typename Tp>
    struct
    log2_sign_t : public log_sign_t<Tp>
    {
      template<typename Up>
	explicit operator Up()
	{ return sign * std::pow(Tp(2), value); };
    };

//  Think of good arithmetic between these types...
