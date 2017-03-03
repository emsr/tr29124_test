/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_ellint test_jacobi_ellint.cpp -lquadmath -Lwrappers/debug -lwrap_boost -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_jacobi_ellint > test_jacobi_ellint.txt

g++ -std=gnu++17 -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_jacobi_ellint test_jacobi_ellint.cpp -lquadmath -Lwrappers/debug -lwrap_boost -lwrap_gsl
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_jacobi_ellint > test_jacobi_ellint.txt
*/

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>

#include "wrap_boost.h"
#include "wrap_gsl.h"

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

void
test_gsl()
{
  std::cout.precision(std::numeric_limits<double>::digits10);
  std::cout << std::showpoint << std::scientific;
  auto width = 8 + std::cout.precision();

  std::vector<double> kvals{0.0, 0.5, 0.75, 0.95, 1.0};

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic sine amplitude function sn(k,u)\n";
  for (auto k : kvals)
    for (int i = -50; i <= 50; ++i)
      {
	auto u = 0.2 * i;
	auto fgnu = __gnu_cxx::jacobi_sn(k, u);
	auto fgsl = gsl::jacobi_sn(k, u);
	std::cout << ' ' << std::setw(width) << u;
	  std::cout << ' ' << std::setw(width) << fgnu
		    << ' ' << std::setw(width) << fgsl
		    << ' ' << fgnu - fgsl;
	std::cout << '\n';
      }

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic cosine amplitude function cn(k,u)\n";
  for (auto k : kvals)
    for (int i = -50; i <= 50; ++i)
      {
	auto u = 0.2 * i;
	auto fgnu = __gnu_cxx::jacobi_cn(k, u);
	auto fgsl = gsl::jacobi_cn(k, u);
	std::cout << ' ' << std::setw(width) << u;
	  std::cout << ' ' << std::setw(width) << fgnu
		    << ' ' << std::setw(width) << fgsl
		    << ' ' << fgnu - fgsl;
	std::cout << '\n';
      }

  std::cout << "\n\n";
  std::cout << "Jacobi elliptic delta amplitude function dn(k,u)\n";
  for (auto k : kvals)
    for (int i = -50; i <= 50; ++i)
      {
	auto u = 0.2 * i;
	auto fgnu = __gnu_cxx::jacobi_dn(k, u);
	auto fgsl = gsl::jacobi_dn(k, u);
	std::cout << ' ' << std::setw(width) << u;
	  std::cout << ' ' << std::setw(width) << fgnu
		    << ' ' << std::setw(width) << fgsl
		    << ' ' << fgnu - fgsl;
	std::cout << '\n';
      }
}

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

  test_gsl();
}
