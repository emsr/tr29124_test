#include <vector>
#include <complex>
#include <iostream>

#include <emsr/traits_utils.h>

int
main()
{
  std::cout << "has_value_type\n";
  std::cout << std::boolalpha;

  auto c = emsr::detail::has_value_type_v<std::complex<double>>;
  std::cout << "std::complex<double>: " << c << '\n';

  auto ci = emsr::detail::has_value_type_v<int>;
  std::cout << "int                 : " << ci << '\n';

  auto ca = emsr::detail::has_value_type_v<double [3]>;
  std::cout << "double a[3]         : " << ca << '\n';

  auto s = std::is_same_v<double, emsr::value_t<std::complex<double>>>;
  std::cout << "double == value_t<std::complex<double>>: " << s << '\n';

  auto si = std::is_same_v<int, emsr::value_t<int>>;
  std::cout << "int    == value_t<int>                 : " << si << '\n';

  auto sa = std::is_same_v<double, emsr::value_t<double [3]>>;
  std::cout << "double == value_t<double [3]>          : " << sa << '\n';
}
