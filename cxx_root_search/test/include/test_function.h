#ifndef ROOT_SEARCH_TEST_FUNCTION_H
#define ROOT_SEARCH_TEST_FUNCTION_H 1

  template<typename Tp>
    constexpr Tp pi = Tp{3.1415926535897932384626433832795028841971L};

  template<typename Tp>
    constexpr Tp e = Tp{2.7182818284590452353602874713526624977572L};

  template<typename Tp>
    constexpr Tp eps = std::sqrt(std::numeric_limits<Tp>::epsilon());

  // f(x) = x^{20} - 1
  // f'(x) = 20 x^{19}
  // f''(x) = 380 x^{18}
  // zero at x = 1 or -1
  template<typename Tp>
    inline Tp
    func1(Tp x)
    { return std::pow(x, Tp{20}) - Tp{1}; }

  template<typename Tp>
    inline Tp
    func1_df(Tp x)
    { return Tp{20} * std::pow(x, Tp{19}); }

  template<typename Tp>
    inline Tp
    func1_df2(Tp x)
    { return Tp{20} * Tp{19} * std::pow(x, Tp{18}); }

  template<typename Tp>
    inline std::pair<Tp, Tp>
    func1_fdf(Tp x)
    { return {func1(x), func1_df(x)}; }

  // f(x) = signum(x) sqrt(std::abs(x))
  // f'(x) = 1 / sqrt(std::abs(x))
  // f''(x) = -signum(x) / pow(abs(x), 3/2) / 2
  // zero at x = 0
  template<typename Tp>
    inline Tp
    func2(Tp x)
    {
      Tp delta;

      if (x > Tp{0})
	delta = Tp{1};
      else if (x < Tp{0})
	delta = Tp{-1};
      else
	delta = Tp{0};

      return std::sqrt(std::abs(x)) * delta;
    }

  template<typename Tp>
    inline Tp
    func2_df(Tp x)
    { return Tp{1} / std::sqrt(std::abs(x)); }

  template<typename Tp>
    inline Tp
    func2_df2(Tp x)
    {
      Tp delta;

      if (x == Tp{0})
	return Tp{0};
      else
	return Tp(x < 0 ? +1 : -1)
		/ std::sqrt(std::abs(x)) / std::abs(x) / Tp{2};
    }

  template<typename Tp>
    inline std::pair<Tp, Tp>
    func2_fdf(Tp x)
    { return {func2(x), func2_df(x)}; }

  // f(x) = x^2 - 1e-8
  // f'(x) = 2x
  // f''(x) = 2
  // zero at x = sqrt(1e-8) or -sqrt(1e-8)
  template<typename Tp>
    inline Tp
    func3(Tp x)
    { return std::pow(x, Tp{2}) - Tp{1e-8}; }

  template<typename Tp>
    inline Tp
    func3_df(Tp x)
    { return Tp{2} * x; }

  template<typename Tp>
    inline Tp
    func3_df2(Tp)
    { return Tp{2}; }

  template<typename Tp>
    inline std::pair<Tp, Tp>
    func3_fdf(Tp x)
    { return {func3(x), func3_df(x)}; }

  // f(x) = x exp(-x)
  // f'(x) = (1 - x) exp(-x)
  // f''(x) = (x - 2) exp(-x)
  // zero at x = 0
  template<typename Tp>
    inline Tp
    func4(Tp x)
    { return x * std::exp(-x); }

  template<typename Tp>
    inline Tp
    func4_df(Tp x)
    { return (Tp{1} - x) * std::exp(-x); }

  template<typename Tp>
    inline Tp
    func4_df2(Tp x)
    { return  (x - Tp{2}) * std::exp(-x); }

  template<typename Tp>
    inline std::pair<Tp, Tp>
    func4_fdf(Tp x)
    { return {func4(x), func4_df(x)}; }

  // f(x) = 1 / (1 + exp(x))
  // f'(x) = -exp(x) / (1 + exp(x))^2
  // f''(x) = -exp(x) / (1 + exp(x))^2 + 2 exp(2x) / (1 + exp(x))^3
  // no roots!
  template<typename Tp>
    inline Tp
    func5(Tp x)
    { return Tp{1} / (Tp{1} + std::exp(x)); }

  template<typename Tp>
    inline Tp
    func5_df(Tp x)
    { return -std::exp(x) / std::pow(Tp{1} + std::exp(x), Tp{2}); }

  template<typename Tp>
    inline Tp
    func5_df2(Tp x)
    { return func5_df(x) * (Tp{1} - Tp{2} * std::exp(x) / (Tp{1} + std::exp(x))); }

  template<typename Tp>
    inline std::pair<Tp, Tp>
    func5_fdf(Tp x)
    { return {func5(x), func5_df(x)}; }

  // f(x) = (x - 1)^7
  // f'(x) = 7 (x - 1)^6
  // f''(x) = 42 (x - 1)^5
  // zero at x = 1
  template<typename Tp>
    inline Tp
    func6(Tp x)
    { return std::pow(x - Tp{1}, Tp{7}); }

  template<typename Tp>
    inline Tp
    func6_df(Tp x)
    { return Tp{7} * std::pow(x - Tp{1}, Tp{6}); }

  template<typename Tp>
    inline Tp
    func6_df2(Tp x)
    { return Tp{42} * std::pow(x - Tp{1}, Tp{5}); }

  template<typename Tp>
    inline std::pair<Tp, Tp>
    func6_fdf(Tp x)
    { return {func6(x), func6_df(x)}; }

  template<typename Tp>
    std::tuple<Tp, Tp, Tp>
    func6_fdf2(Tp x)
    { return {func6(x), func6_df(x), func6_df(x)}; }

  // f(x) = sin(x)
  // f'(x) = cos(x)
  // f''(x) = -sin(x)
  template<typename Tp>
    inline Tp
    sin_f(Tp x)
    { return std::sin(x); }

  template<typename Tp>
    inline Tp
    sin_df(Tp x)
    { return std::cos(x); }

  template<typename Tp>
    inline Tp
    sin_df2(Tp x)
    { return -std::sin(x); }

  template<typename Tp>
    inline std::pair<Tp, Tp>
    sin_fdf(Tp x)
    { return {std::sin(x), std::cos(x)}; }

  // f(x) = cos(x)
  // f'(x) = -sin(x)
  // f''(x) = -cos(x)
  template<typename Tp>
    inline Tp
    cos_f(Tp x)
    { return std::cos(x); }

  template<typename Tp>
    inline Tp
    cos_df(Tp x)
    { return -std::sin(x); }

  template<typename Tp>
    inline Tp
    cos_df2(Tp x)
    { return -std::cos(x); }

  template<typename Tp>
    inline std::pair<Tp, Tp>
    cos_fdf(Tp x)
    { return {std::cos(x), -std::sin(x)}; }

  // Linear function to test that solvers exit correctly 
  // when entered with an exact root
  // f(x) = -pi x + e
  // f'(x) = -pi
  // f''(x) = 0
  template<typename Tp>
    inline Tp
    func7(Tp x)
    { return -pi<Tp> * x + e<Tp>; }

  template<typename Tp>
    inline Tp
    func7_df(Tp)
    { return -pi<Tp>; }

  template<typename Tp>
    inline Tp
    func7_df2(Tp)
    { return Tp{0}; }

  template<typename Tp>
    inline std::pair<Tp, Tp>
    func7_fdf(Tp x)
    { return {func7(x), func7_df(x)}; }

  template<typename _Tp>
    using func_value_t = _Tp (*)(_Tp);

  template<typename Tp>
    struct value_test_t
    {
      func_value_t<Tp> func;
      Tp lower;
      Tp upper;
      Tp eps;
      Tp root;
      const char* str;
    };

  // Array of test functions.
  template<typename Tp>
    value_test_t<Tp>
    test_value_funcs[10]
    {
      {sin_f<Tp>, Tp{3}, Tp{4}, eps<Tp>, pi<Tp>, "sin(x)"},
      {sin_f<Tp>, Tp{-4}, Tp{-3}, eps<Tp>, -pi<Tp>, "sin(x)"},
      {sin_f<Tp>, Tp{-1} / Tp{3}, Tp{1}, eps<Tp>, Tp{0}, "sin(x)"},
      {cos_f<Tp>, Tp{0}, Tp{3}, eps<Tp>, pi<Tp> / Tp{2}, "cos(x)"},
      {cos_f<Tp>, Tp{-3}, Tp{0}, eps<Tp>, -pi<Tp> / Tp{2}, "cos(x)"},
      {func1<Tp>, Tp{0.1L}, Tp{2}, eps<Tp>, Tp{1}, "x^20 - 1"},
      {func2<Tp>, Tp{-1} / Tp{3}, Tp{1}, eps<Tp>, Tp{0}, "1/sqrt(abs(x))"},
      {func3<Tp>, Tp{0}, Tp{1}, eps<Tp>, Tp(std::sqrt(1e-8L)), "x^2 - 1e-8"},
      {func4<Tp>, Tp{-1} / Tp{3}, Tp{2}, eps<Tp>, Tp{0}, "x exp(-x)"},
      {func6<Tp>, Tp{0.9995L}, Tp{1.0002L}, eps<Tp>, Tp{1}, "(x - 1)^7"}
    };

  template<typename Tp>
    using func_fdf_t = std::pair<Tp, Tp> (*)(Tp);

  template<typename Tp>
    struct fdf_test_t
    {
      func_fdf_t<Tp> func;
      Tp lower;
      Tp upper;
      Tp trial;
      Tp eps;
      Tp root;
      const char* str;
    };

  // Array of test functions returning {value, deriv}.
  template<typename Tp>
    fdf_test_t<Tp>
    test_fdf_funcs[12]
    {
      {sin_fdf<Tp>, Tp{3}, Tp{4}, Tp{3.4L}, eps<Tp>, pi<Tp>, "sin(x)"},
      {sin_fdf<Tp>, Tp{-4}, Tp{-3}, Tp{-3.3L}, eps<Tp>, -pi<Tp>, "sin(x)"},
      {sin_fdf<Tp>, Tp{-1} / Tp{3}, Tp{1}, Tp{0.5L}, eps<Tp>, Tp{0}, "sin(x)"},
      {cos_fdf<Tp>, Tp{0}, Tp{3}, Tp{0.6L}, eps<Tp>, pi<Tp> / Tp{2}, "cos(x)"},
      {cos_fdf<Tp>, Tp{-3}, Tp{0}, Tp{-2.5L}, eps<Tp>, -pi<Tp> / Tp{2}, "cos(x)"},
      {func1_fdf<Tp>, Tp{0.1L}, Tp{2}, Tp{0.9L}, eps<Tp>, Tp{1}, "x^20 - 1"},
      {func1_fdf<Tp>, Tp{0.1L}, Tp{2}, Tp{1.1L}, eps<Tp>, Tp{1}, "x^20 - 1"},
      {func2_fdf<Tp>, Tp{-1} / Tp{3}, Tp{1}, Tp{0.001}, eps<Tp>, Tp{0}, "1/sqrt(abs(x))"},
      {func3_fdf<Tp>, Tp{0}, Tp{1}, Tp{1}, eps<Tp>, Tp(std::sqrt(1e-8L)), "x^2 - 1e-8"},
      {func4_fdf<Tp>, Tp{-1} / Tp{3}, Tp{2}, Tp{2}, eps<Tp>, Tp{0}, "x exp(-x)"},
      {func5_fdf<Tp>, Tp{-1}, Tp{0.1L}, Tp{0}, eps<Tp>, Tp{0}, "1 / (1 + exp(x))"},
      {func7_fdf<Tp>, Tp{1}, Tp{2}, Tp{1.5L}, eps<Tp>, e<Tp> / pi<Tp>, "-pi x + e"}
    };

  template<typename Tp>
    using func_fdf2_t = std::tuple<Tp, Tp, Tp> (*)(Tp);

  template<typename Tp>
    struct fdf2_test_t
    {
      func_fdf2_t<Tp> func;
      Tp lower;
      Tp upper;
      Tp trial;
      Tp eps;
      Tp root;
      const char* str;
    };

  // Array of test functions returning {value, deriv, deriv2}.
  template<typename Tp>
    fdf2_test_t<Tp>
    test_fdf2_funcs[1]
    {
      {func6_df2<Tp>, Tp{0.9995L}, Tp{1.0002L}, Tp{0.9L}, eps<Tp>, Tp{1}, "(x - 1)^7"}
    };

#endif // ROOT_SEARCH_TEST_FUNCTION_H
