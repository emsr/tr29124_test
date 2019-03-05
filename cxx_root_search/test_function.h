#ifndef CXX_ROOT_SEARCH_TEST_FUNCTION_H
#define CXX_ROOT_SEARCH_TEST_FUNCTION_H 1

  // f(x) = x^{20} - 1
  // f'(x) = 20x^{19}
  // zero at x = 1 or -1
  template<typename _Tp>
    _Tp func1(_Tp x);

  template<typename _Tp>
    _Tp func1_df(_Tp x);

  template<typename _Tp>
    std::pair<_Tp, _Tp> func1_fdf(_Tp x);

  // f(x) = std::sqrt(std::abs(x))*sgn(x)
  // f'(x) = 1 / std::sqrt(std::abs(x)
  // zero at x = 0
  template<typename _Tp>
    _Tp func2(_Tp x);

  template<typename _Tp>
    _Tp func2_df(_Tp x);

  template<typename _Tp>
    std::pair<_Tp, _Tp> func2_fdf(_Tp x);

  // f(x) = x^2 - 1e-8
  // f'(x) = 2x
  // zero at x = std::sqrt(1e-8) or -std::sqrt(1e-8)
  template<typename _Tp>
    _Tp func3(_Tp x);

  template<typename _Tp>
    _Tp func3_df(_Tp x);

  template<typename _Tp>
    std::pair<_Tp, _Tp> func3_fdf(_Tp x);

  // f(x) = x std::exp(-x)
  // f'(x) = std::exp(-x) - x std::exp(-x)
  // zero at x = 0
  template<typename _Tp>
    _Tp func4(_Tp x);

  template<typename _Tp>
    _Tp func4_df(_Tp x);

  template<typename _Tp>
    std::pair<_Tp, _Tp> func4_fdf(_Tp x);

  // f(x) = 1/(1+std::exp(x))
  // f'(x) = -std::exp(x) /(1 + std::exp(x))^2
  // no roots!
  template<typename _Tp>
    _Tp func5(_Tp x);

  template<typename _Tp>
    _Tp func5_df(_Tp x);

  template<typename _Tp>
    std::pair<_Tp, _Tp> func5_fdf(_Tp x);

  // f(x) = (x - 1)^7
  // f'(x) = 7 * (x - 1)^6
  // zero at x = 1
  template<typename _Tp>
    _Tp func6(_Tp x);

  template<typename _Tp>
    _Tp func6_df(_Tp x);

  template<typename _Tp>
    std::pair<_Tp, _Tp> func6_fdf(_Tp x);

  template<typename _Tp>
    _Tp func6_df2(_Tp x);

  template<typename _Tp>
    std::tuple<_Tp, _Tp, _Tp> func6_fdf2(_Tp x);

  // f(x) = sin(x)
  template<typename _Tp>
    _Tp sin_f(_Tp x);

  template<typename _Tp>
    _Tp sin_df(_Tp x);

  template<typename _Tp>
    std::pair<_Tp, _Tp> sin_fdf(_Tp x);

  // f(x) = cos(x)
  template<typename _Tp>
    _Tp cos_f(_Tp x);

  template<typename _Tp>
    _Tp cos_df(_Tp x);

  template<typename _Tp>
    std::pair<_Tp, _Tp> cos_fdf(_Tp x);

  // linear function to test that solvers exit correctly 
  // when entered with an exact root
  template<typename _Tp>
    _Tp func7(_Tp x);

  template<typename _Tp>
    _Tp func7_df(_Tp x);

  template<typename _Tp>
    std::pair<_Tp, _Tp> func7_fdf(_Tp x);

  template<typename _Tp>
    using func_value_t = _Tp (*)(_Tp);

  template<typename _Tp>
    struct value_test_t
    {
      func_value_t<_Tp> func;
      _Tp lower;
      _Tp upper;
      _Tp eps;
      _Tp root;
      const char* str;
    };

  template<typename _Tp>
    constexpr _Tp pi = _Tp{3.1415926535897932384626433832795028841971L};

  template<typename _Tp>
    constexpr _Tp e = _Tp{2.7182818284590452353602874713526624977572L};

  template<typename _Tp>
    constexpr _Tp eps = std::sqrt(std::numeric_limits<_Tp>::epsilon());

  // Array of test functions.
  template<typename _Tp>
    value_test_t<_Tp>
    test_value_funcs[10]
    {
      {sin_f<_Tp>, _Tp{3}, _Tp{4}, eps<_Tp>, pi<_Tp>, "sin(x)"},
      {sin_f<_Tp>, _Tp{-4}, _Tp{-3}, eps<_Tp>, -pi<_Tp>, "sin(x)"},
      {sin_f<_Tp>, _Tp{-1} / _Tp{3}, _Tp{1}, eps<_Tp>, _Tp{0}, "sin(x)"},
      {cos_f<_Tp>, _Tp{0}, _Tp{3}, eps<_Tp>, pi<_Tp> / _Tp{2}, "cos(x)"},
      {cos_f<_Tp>, _Tp{-3}, _Tp{0}, eps<_Tp>, -pi<_Tp> / _Tp{2}, "cos(x)"},
      {func1<_Tp>, _Tp{0.1L}, _Tp{2}, eps<_Tp>, _Tp{1}, "x^20 - 1"},
      {func2<_Tp>, _Tp{-1} / _Tp{3}, _Tp{1}, eps<_Tp>, _Tp{0}, "1/sqrt(abs(x))"},
      {func3<_Tp>, _Tp{0}, _Tp{1}, eps<_Tp>, _Tp(std::sqrt(1e-8L)), "x^2 - 1e-8"},
      {func4<_Tp>, _Tp{-1} / _Tp{3}, _Tp{2}, eps<_Tp>, _Tp{0}, "x exp(-x)"},
      {func6<_Tp>, _Tp{0.9995L}, _Tp{1.0002L}, eps<_Tp>, _Tp{1}, "(x - 1)^7"}
    };

  template<typename _Tp>
    using func_fdf_t = std::pair<_Tp, _Tp> (*)(_Tp);

  template<typename _Tp>
    struct fdf_test_t
    {
      func_fdf_t<_Tp> func;
      _Tp lower;
      _Tp upper;
      _Tp trial;
      _Tp eps;
      _Tp root;
      const char* str;
    };

  // Array of test functions returning {value, deriv}.
  template<typename _Tp>
    fdf_test_t<_Tp>
    test_fdf_funcs[12]
    {
      {sin_fdf<_Tp>, _Tp{3}, _Tp{4}, _Tp{3.4L}, eps<_Tp>, pi<_Tp>, "sin(x)"},
      {sin_fdf<_Tp>, _Tp{-4}, _Tp{-3}, _Tp{-3.3L}, eps<_Tp>, -pi<_Tp>, "sin(x)"},
      {sin_fdf<_Tp>, _Tp{-1} / _Tp{3}, _Tp{1}, _Tp{0.5L}, eps<_Tp>, _Tp{0}, "sin(x)"},
      {cos_fdf<_Tp>, _Tp{0}, _Tp{3}, _Tp{0.6L}, eps<_Tp>, pi<_Tp> / _Tp{2}, "cos(x)"},
      {cos_fdf<_Tp>, _Tp{-3}, _Tp{0}, _Tp{-2.5L}, eps<_Tp>, -pi<_Tp> / _Tp{2}, "cos(x)"},
      {func1_fdf<_Tp>, _Tp{0.1L}, _Tp{2}, _Tp{0.9L}, eps<_Tp>, _Tp{1}, "x^20 - 1"},
      {func1_fdf<_Tp>, _Tp{0.1L}, _Tp{2}, _Tp{1.1L}, eps<_Tp>, _Tp{1}, "x^20 - 1"},
      {func2_fdf<_Tp>, _Tp{-1} / _Tp{3}, _Tp{1}, _Tp{0.001}, eps<_Tp>, _Tp{0}, "1/sqrt(abs(x))"},
      {func3_fdf<_Tp>, _Tp{0}, _Tp{1}, _Tp{1}, eps<_Tp>, _Tp(std::sqrt(1e-8L)), "x^2 - 1e-8"},
      {func4_fdf<_Tp>, _Tp{-1} / _Tp{3}, _Tp{2}, _Tp{2}, eps<_Tp>, _Tp{0}, "x exp(-x)"},
      {func5_fdf<_Tp>, _Tp{-1}, _Tp{0.1L}, _Tp{0}, eps<_Tp>, _Tp{0}, "1 / (1 + exp(x))"},
      {func7_fdf<_Tp>, _Tp{1}, _Tp{2}, _Tp{1.5L}, eps<_Tp>, e<_Tp> / pi<_Tp>, "-pi x + e"}
    };

  template<typename _Tp>
    using func_fdf2_t = std::tuple<_Tp, _Tp, _Tp> (*)(_Tp);

  template<typename _Tp>
    struct fdf2_test_t
    {
      func_fdf2_t<_Tp> func;
      _Tp lower;
      _Tp upper;
      _Tp trial;
      _Tp eps;
      _Tp root;
      const char* str;
    };

  // Array of test functions returning {value, deriv, deriv2}.
  template<typename _Tp>
    fdf2_test_t<_Tp>
    test_fdf2_funcs[1]
    {
      {func6_df2<_Tp>, _Tp{0.9995L}, _Tp{1.0002L}, _Tp{0.9L}, eps<_Tp>, _Tp{1}, "(x - 1)^7"}
    };

#endif // CXX_ROOT_SEARCH_TEST_FUNCTION_H

#include "test_function.tcc"
