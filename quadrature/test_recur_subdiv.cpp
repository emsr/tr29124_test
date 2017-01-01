/*
$HOME/bin/bin/g++-std=c++17 -Wall -Wextra -I.. -o test_recur_subdiv test_recur_subdiv.cpp -lquadmath
./test_recur_subdiv
*/

#include <iostream>
#include <iomanip>

#include <experimental/array>
#include <vector>
#include <cmath>

template<typename _Tp>
  struct __subdiv
  {
    std::size_t __n_segs;
    std::vector<_Tp> __pt;
    std::vector<std::array<std::size_t, 2>> __seg;

    __subdiv(std::size_t __n_pts, _Tp __a, _Tp __b);

    __subdiv(std::size_t __n_pts, _Tp __a, _Tp __b, std::size_t n_levels);

    std::size_t
    n_levels() const
    { return std::log2(__seg.size() / __n_segs); }

  };

/**
 * 
 */
template<typename _Tp>
  __subdiv<_Tp>::__subdiv(std::size_t __n_pts, _Tp __a, _Tp __b)
  : __n_segs(__n_pts - 1),
    __pt(__n_pts),
    __seg(__n_segs)
  {
    auto __delta = (__b - __a) / __n_pts - 1;
    for (auto __i = 0u; __i < __n_segs; ++__i)
      __pt[__i] = __a + __i * __delta;
    __pt[__n_segs] = __b;

    for (auto __i = 0u; __i < __n_segs; ++__i)
      {
	__seg[__i][0] = __i;
	__seg[__i][1] = __i + 1;
      }
  }

/**
 * 
 */
template<typename _Tp>
  __subdiv<_Tp>::__subdiv(std::size_t __n_pts, _Tp __a, _Tp __b,
			  std::size_t __n_levels)
  : __n_segs{__n_pts - 1},
    __pt{},
    __seg{}
  {
    std::size_t __size = __n_segs;
    for (std::size_t __l = 0u; __l < __n_levels; ++__l)
      __size *= 2;

    __pt.reserve(__size + 1);
    __seg.reserve(__size);
    auto __delta = (__b - __a) / __n_pts - 1;
    for (std::size_t __i = 0u; __i < __n_segs; ++__i)
      __pt.push_back(__a + __i * __delta);
    __pt.push_back(__b);

    for (std::size_t __i = 0u; __i < __n_segs; ++__i)
      __seg.push_back(std::experimental::make_array(__i, __i + 1));

    auto __start = 0u;
    auto __stop = __n_segs;
    for (std::size_t __l = 0u; __l < __n_levels; ++__l)
      {
	for (std::size_t __i = __start; __i < __stop; ++__i)
	  {
	    auto __s0 = __seg[__i][0];
	    auto __s1 = __seg[__i][1];
	    __pt.push_back((__pt[__s0] + __pt[__s1]) / _Tp{2});
	    auto __p = ++__n_pts;
	    __seg.push_back(std::experimental::make_array(__s0, __p));
	    __seg.push_back(std::experimental::make_array(__p, __s1));
	  }
	__start = __stop;
	__stop *= 2;
      }
  }

template<typename _Tp>
  void
  test_recur_subdiv(std::size_t n_pts, _Tp a, _Tp b, std::size_t n_levels)
  {
    __subdiv<_Tp> sdiv(n_pts, a, b, n_levels);
    std::cout << "\npoints\n";
    for (auto pt : sdiv.__pt)
      std::cout << ' ' << std::setw(6) << pt << '\n';
    std::cout << "\nsegments\n";
    for (const auto& seg : sdiv.__seg)
      std::cout << ' ' << std::setw(4) << seg[0]
	       << ", " << std::setw(4) << seg[1] << '\n';
  }

int
main()
{
  test_recur_subdiv(10, -1.0, +1.0, 5);
}
