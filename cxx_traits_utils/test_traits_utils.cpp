#include <vector>
#include <complex>
#include <iostream>

#include <ext/traits_utils.h>

int
main()
{
  std::cout << "has_value_type\n";
  std::cout << std::boolalpha;

  auto c = __gnu_cxx::__detail::__has_value_type_v<std::complex<double>>;
  std::cout << "std::complex<double>: " << c << '\n';

  auto ci = __gnu_cxx::__detail::__has_value_type_v<int>;
  std::cout << "int                 : " << ci << '\n';

  auto ca = __gnu_cxx::__detail::__has_value_type_v<double [3]>;
  std::cout << "double a[3]         : " << ca << '\n';

  auto s = std::is_same_v<double, __gnu_cxx::__value_t<std::complex<double>>>;
  std::cout << "double == value_t<std::complex<double>>: " << s << '\n';

  auto si = std::is_same_v<int, __gnu_cxx::__value_t<int>>;
  std::cout << "int    == value_t<int>                 : " << si << '\n';

  auto sa = std::is_same_v<double, __gnu_cxx::__value_t<double [3]>>;
  std::cout << "double == value_t<double [3]>          : " << sa << '\n';
}
