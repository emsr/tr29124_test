
#include <complex>
#include <complex_algorithm.h>

int
main()
{
#if __cpp_lib_interpolate >=  201902L
  std::complex<double> a{1, 2}, b{3, 4};
  auto l = std::lerp(a, b, 0.3);
#endif // __cpp_lib_interpolate >=  201902L
}
