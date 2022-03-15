
#include <complex>
#include <emsr/complex_util.h>

template<typename _Tp>
  struct
  test_is_complex_neg
  {
    static_assert(!emsr::is_complex_v<_Tp>);
    static_assert(!emsr::is_complex_v<const _Tp>);
    static_assert(!emsr::is_complex_v<volatile _Tp>);
  };

template<typename _Tp>
  struct
  test_is_complex
  {
    static_assert(emsr::is_complex_v<_Tp>);
    static_assert(emsr::is_complex_v<const _Tp>);
    static_assert(emsr::is_complex_v<volatile _Tp>);
  };

test_is_complex_neg<int> i;

test_is_complex_neg<float> f;
test_is_complex_neg<double> d;
test_is_complex_neg<long double> ld;

test_is_complex<std::complex<float>> cf;
test_is_complex<std::complex<double>> cd;
test_is_complex<std::complex<long double>> cld;

int
main()
{
}
