/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_ellint test_jacobi_ellint.cpp -L$HOME/bin/lib64 -lquadmath
./test_jacobi_ellint > test_jacobi_ellint.txt

g++ -std=gnu++17 -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_ellint test_jacobi_ellint.cpp -lquadmath
./test_jacobi_ellint > test_jacobi_ellint.txt
*/

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_sc(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::__detail::__jacobi_sncndn<__type>(__k, __u).sc();
    }

  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_sd(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::__detail::__jacobi_sncndn<__type>(__k, __u).sd();
    }

  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_cd(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::__detail::__jacobi_sncndn<__type>(__k, __u).cd();
    }

  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_cs(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::__detail::__jacobi_sncndn<__type>(__k, __u).cs();
    }

  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_dc(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::__detail::__jacobi_sncndn<__type>(__k, __u).dc();
    }

  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_ds(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::__detail::__jacobi_sncndn<__type>(__k, __u).ds();
    }

  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_nc(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::__detail::__jacobi_sncndn<__type>(__k, __u).nc();
    }

  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_nd(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::__detail::__jacobi_sncndn<__type>(__k, __u).nd();
    }

  template<typename _Kp, typename _Up>
    inline __gnu_cxx::__promote_fp_t<_Kp, _Up>
    jacobi_ns(_Kp __k, _Up __u)
    {
      using __type = __gnu_cxx::__promote_fp_t<_Kp, _Up>;
      return std::__detail::__jacobi_sncndn<__type>(__k, __u).ns();
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

int
main()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();

  std::vector<double> kvals{0.0, 0.5, 0.75, 0.95, 1.0};

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic sine amplitude function sn(k,u)\n";
  for (int i = -200; i <= 200; ++i)
    {
      auto u = 0.05 * i;
      std::cout << ' ' << std::setw(width) << u;
      for (auto k : kvals)
	std::cout << ' ' << std::setw(width) << __gnu_cxx::jacobi_sn(k, u);
      std::cout << '\n';
    }

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic cosine amplitude function cn(k,u)\n";
  for (int i = -200; i <= 200; ++i)
    {
      auto u = 0.05 * i;
      std::cout << ' ' << std::setw(width) << u;
      for (auto k : kvals)
	std::cout << ' ' << std::setw(width) << __gnu_cxx::jacobi_cn(k, u);
      std::cout << '\n';
    }

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic delta amplitude function dn(k,u)\n";
  for (int i = -200; i <= 200; ++i)
    {
      auto u = 0.05 * i;
      std::cout << ' ' << std::setw(width) << u;
      for (auto k : kvals)
	std::cout << ' ' << std::setw(width) << __gnu_cxx::jacobi_dn(k, u);
      std::cout << '\n';
    }
}
