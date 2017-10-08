/*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_polynomial_root test_polynomial_root.cpp -lquadmath
./test_polynomial_root > test_polynomial_root.txt
*/

#include <complex>
#include <vector>
#include <iostream>

#include "polynomial.h"

template<typename _Tp>
  int
  polynomial_roots(const __gnu_cxx::_Polynomial<_Tp>& poly,
		   std::vector<std::complex<_Tp>>& root)
  {
    auto xcof = poly;
    __gnu_cxx::_Polynomial<_Tp> cof(xcof);

    _Tp mag, cofj;
    std::complex<_Tp> t, t1, u, ud;

    bool final = false;
    auto n = xcof.degree() - 1;
    if (n <= 0)
      return 1;
    if (n > 36)
      return 2;
    if (xcof.coefficient(n) == _Tp{0})
      return 4;

    //int retry = 0;
    auto n1 = n;
    auto n2 = n;
    auto nsav = n;
    for (std::size_t j = 0; j <= nsav; ++j)
      cof.coefficient(n - j, xcof.coefficient(j));

    std::complex<_Tp> xsav;

    while (n > 0)
      {
	std::complex<_Tp> x0{0.00500101, 0.01000101};
	int retry = 0;

  TRY_AGAIN:

	++retry;
	// It would be nice to have a complex tie...
	//auto& x0r = x0.real();
	//auto& x0i = x0.imag();
	//x0i = _Tp{-10} * std::exchange(x0r, _Tp{-10} * x0i);
	auto x0i = x0.imag();
	x0.imag(_Tp{-10} * x0.real());
	x0.real(_Tp{-10} * x0i);
	auto x = x0;

  FINAL_ITER:

	int iter = 0;

	bool real_root = false;
	while (iter < 500)
	{
	  u.real(cof.coefficient(n));
	  if (u.real() == _Tp{0})
	    { // This root is zero
	      x.real(_Tp{0});
	      --n1;
	      --n2;
	      goto ZERO_ROOT;
	    }
	  u.imag(_Tp{0});
	  auto ud = std::complex<_Tp>{};
	  auto t = std::complex<_Tp>{1, 0};
	  for (std::size_t i = 0; i < n; ++i) 
	    {
	      auto t1 = x * t;
	      cofj = cof.coefficient(n - 1 - i); // Evaluate polynomial
	      u += cofj * t1;
	      cofj *= _Tp(i + 1); // Evaluate derivative
	      ud += cofj * conj(t);
	      t = t1;
	    }

	  auto mag = std::norm(ud);
	  if (mag == _Tp{0}) 
	    {
	      if (!final)
		goto TRY_AGAIN;
	      x = xsav;
	      goto findon;
	    }
	  auto dx = u * ud / mag;
	  x += dx;
	  if (std::abs(dx.imag()) + std::abs(dx.real()) < 1.0e-6) 
	    goto LOOP_DONE;
	  ++iter;
	} // while iter < 500

	if (final)
	  goto LOOP_DONE;
	if (retry < 5)
	  goto TRY_AGAIN;

	return 3;

  LOOP_DONE:

	// Swap original and reduced polynomials
	for (std::size_t j = 0; j <= n2; ++j)
	  {
	    auto poo = cof.coefficient(j);
	    cof.coefficient(j, xcof.coefficient(nsav - j));
	    xcof.coefficient(nsav - j, poo);
	  }
	n = std::exchange(n1, n);

	if (!final)
	  {
	    final = true;
	    if (std::abs(x.imag() / x.real()) < 1.0e-4) 
	      x.imag(_Tp{0});
	    xsav = x;
	    goto FINAL_ITER; // do final iteration on original polynomial
	  }

  findon:

	final = false;
	if (std::abs(x.imag() / x.real()) >= 1.0e-5)
	  {
	    cofj = _Tp{2} * x.real();
	    mag = std::norm(x);
	    n -= 2;
	  }
	else
	  { // root is real
  ZERO_ROOT:
	    real_root = true;
	    x.imag(_Tp{0});
	    cofj = x.real();
	    mag = _Tp{0};
	    --n;
	  }

	// Divide working polynomial cof(z) by z - x
	cof.coefficient(1, cof.coefficient(1) + cofj * cof.coefficient(0));
	for (std::size_t j = 1; j < n; ++j) 
	  cof.coefficient(j + 1, cof.coefficient(j + 1) + cofj * cof.coefficient(j) - mag * cof.coefficient(j - 1));

	root.push_back(x);
	if (!real_root) 
	  root.push_back(std::conj(x));
      }

    return 0;
  }


int
main()
{
  __gnu_cxx::_Polynomial<double> P({1.0, -2.0, 1.0, -4.0, 2,0});
  std::vector<std::complex<double>> zero;
  polynomial_roots(P, zero);
  std::cout << "Found " << zero.size() << " roots:\n";
  for (std::size_t i = 0; i < zero.size(); ++i)
    std::cout << "Root " << (i + 1) << ": " << zero[i] << '\n';
}

