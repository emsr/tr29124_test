/*
g++ -Wall -Wextra -o oblate_harmonic oblate_harmonic.cpp
./oblate_harmonic > oblate_harmonic.txt
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

double factco(int i, double pn, int m, double huge);
double factconew(int n, double pn, int m);

/**
 *  Calculation of oblate spheroidal harmonics
 *
 *    x        Argument of the functions
 *    m        Degree of the spheroidal harmonics
 *    nmax     Maximum order of the functions:
 *   	       We shall get  functions of all the orders below
 *   	       min(nmax,nuevo). Nuevo is defined below.
 *    mode     Mode of calculation (see below)
 *
 * If mode is equal to 0:
 *    Rl[l+1]
 *    Tl[l+1]
 *    nuevo    Maximum order of functions calculated when
 *             Rl[nmax+1] is larger than 1/tiny.
 *             Tiny is defined below.
 *
 * If mode is equal to 1:
 *    Rl[l+1]/(2m-1)!!
 *    Tl[l+1]/(2m)!!
 *    nuevo    Maximum order of functions calculated when
 *             Rl[nmax+1]/(2m-1)!! is larger than 1/tiny.
 *             Tiny is defined below.
 *
 * If mode is equal to 2:
 *    Rl[l+1]*tiny/(2m-1)!!/(x^2+1)^(m/2)
 *    Tl[l+1]*(x^2+1)^(m/2)/(2m)!!/tiny
 *    nuevo    Maximum order of functions calculated when
 *             Rl[nmax+1]/(2m-1)!! is larger than 1/tiny.
 *             Tiny is defined below.
 */
void
oblate_harmonic(double x, int m, int nmax, int mode, std::vector<double> &Rl, std::vector<double> &Tl, int &nuevo)
{
    // dimension Rl(0:nmax+1),Tl(0:nmax+1);
    // Required accuracy for the continued fraction (Lentz-Thompson)
    constexpr double eps = 1.0e-16;
    // Small parameter close to the underflow limit
    constexpr double tiny = 1.0e-290;
    constexpr auto huge = 1.0 / tiny;
    const auto sqrttiny = std::sqrt(tiny);
    if (x <= 0.0)
    {
        std::cout << "improper argument. x must be greater than 0";
        return;
    }
    const auto fl = m / 2.0;
    double ar;
    if (double(int(fl)) == fl)
        ar = 1.0;
    else
        ar = -1.0;
    const auto dzm = std::sqrt(x * x + 1.0);
    // We use the code if nmax is greater than or equal to 2
    auto nmaxp = nmax;
    if (nmax < 2)
        nmax = 2;
    // We store the values of Rl(0) and Rl(1)
    auto rl0 = 1.0;
    if (mode == 0)
    {
        if (m > 0)
        {
            int j = 1;
            for (int i = 1; i <= m; ++i)
            {
                rl0 *= dzm * (2 * m - j);
                j += 2;
            }
        }
        Rl[0] = rl0;
        Rl[1] = x * (2 * m + 1) * Rl[0];
    }
    else if (mode == 1)
    {
        if (m * std::log(dzm) > std::log(huge))
        {
            std::cout << "better try mode=2";
            return;
        }
        if (m > 0)
        {
            for (int i = 1; i <= m; ++i)
            {
                rl0 *= dzm;
            }
        }
        Rl[0] = rl0;
        Rl[1] = x * (2 * m + 1) * Rl[0];
    }
    else
    {
        Rl[0] = rl0 * tiny;
        Rl[1] = x * (2 * m + 1) * Rl[0];
    }
    // We use the recurrence relations
    int np = 1;
    while (np <= nmax && std::abs((np + 1.0) * Rl[np]) < huge)
    {
        Rl[np + 1] = ((2.0 * (np + m) + 1.0) * x * Rl[np] + (np + m + m) * Rl[np - 1]) / (np + 1.0);
        ++np;
    }
    nmax = np - 1;
    double factor;
    if (mode == 0)
    {
        factor = factco(nmax, Rl[nmax + 1], m, huge);
        while (factor == 0.0 && nmax > 0)
        {
            nmax = nmax - 3;
            factor = factco(nmax, Rl[nmax + 1], m, huge);
        }
        if (nmax <= 0)
        {
            std::cout << "try another m";
            return;
        }
    }
    else
    {
        factor = factconew(nmax, Rl[nmax + 1], m);
    }
    nuevo = nmax;
    // We evaluate the continued fraction using Lentz-Thompson algorithm
    auto n = nmax;
    auto mm = 0;
    auto b = (1.0 + double(n + 2) / double(2 * m + n + 1)) * x;
    auto a = 1.0;
    auto fc = sqrttiny;
    auto c0 = fc;
    auto d0 = 0.0;
_10:
    d0 = b + a * d0;
    if (std::abs(d0) < sqrttiny)
        d0 = sqrttiny;
    c0 = b + a / c0;
    if (std::abs(c0) < sqrttiny)
        c0 = sqrttiny;
    d0 = 1.0 / d0;
    auto delta = c0 * d0;
    fc *= delta;
    ++mm;
    a = double(n + mm + 1) / double(n + 2 * m + mm);
    b = (1.0 + double(n + mm + 2) / double(2 * m + n + mm + 1)) * x;
    if (std::abs(delta - 1.0) > eps)
        goto _10;
    // We evaluate Tl[nmax+1] and Tl[nmax] using :
    // The Wronskian : W{Rl(nmax),Tl[nmax]} = 1 / (1 - x^2)
    // The known values of Rl(nmax+1) and Rl(nmax)
    // The value of h = Tl[nmax+1]/Tl[nmax]
    Tl[nmax] = ar * factor / (1.0 + fc * Rl[nmax] / Rl[nmax + 1]);
    Tl[nmax + 1] = Tl[nmax] * fc;
    for (int i = 1; i <= nmax; ++i)
    {
        np = nmax + 1 - i;
        n = np - 1;
        Tl[n] = ((2 * np + 2 * m + 1) * x * Tl[np] + (np + 1) * Tl[np + 1]) / double(np + m + m);
    }
    nmax = nmaxp;
    return;
}

double
factco(int i, double pn, int m, double huge)
{
    double fact = 1.0 / pn;
    if (m > 0)
    {
	int j = 2 * m;
	while (j > 1 && fact < huge)
        {
            fact *= double(i + j);
            --j;
	}
	if (j > 2)
            fact = 0.0;
    }
    else
    {
	fact = 1.0 / (i + 1.0) / pn;
    }
    return fact;
}

double
factconew(int n, double pn, int m)
{
    double fact = 1.0 / (n + 1.0) / pn;
    if (m > 0)
    {
	int j = 2 * m;
	for (int l = 1; l <= n; ++l)
            fact *= double(j + l) / double(l);
    }
    return fact;
}

int
main()
{
  unsigned int m = 3;
  unsigned int n_max = 10;
  int mode = 0, nuevo = 0;
  std::vector<double> R(n_max + 2);
  std::vector<double> T(n_max + 2);
  for (double x : {0.1, 0.5, 1.0, 2.0, 5.0})
  {
    oblate_harmonic(x, m, n_max, mode, R, T, nuevo);
    std::cout << "\nx = " << x;
    std::cout << '\n';
    for (unsigned int n = 0; n <= n_max; ++n)
    {
      std::cout << ' ' << std::setw(12) << R[n];
    }
    std::cout << '\n';
    for (unsigned int n = 0; n <= n_max; ++n)
    {
      std::cout << ' ' << std::setw(12) << T[n];
    }
  }
  std::cout << '\n';
}
