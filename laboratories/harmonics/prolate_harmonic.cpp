/*
g++ -Wall -Wextra -o prolate_harmonic prolate_harmonic.cpp
./prolate_harmonic > prolate_harmonic.txt
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

double factco(int i, double pn, int m, double huge);
double factconew(int n, double pn, int m);

/**
 * Calculation of prolate spheroidal harmonics
 *
 *    x        Argument of the functions
 *    m        Degree of the spheroidal harmonics
 *    nmax     Maximum order of the functions:
 *   	       We shall get functions of all the orders below
 *   	       min(nmax,nuevo). Nuevo is defined below.
 *    mode     Mode of calculation (see below)
 *
 * If mode is equal to 0:
 *    Pl[l+1]
 *    Ql[l+1]
 *    nuevo    Maximum order of  functions calculated when
 *             Pl[nmax+1] is larger than 1/tiny.
 *             Tiny is defined below.
 *
 * If mode is equal to 1:
 *    Pl[l+1]/(2m-1)!!
 *    Ql[l+1]/(2m)!!
 *    nuevo    Maximum order of functions calculated when
 *             Pl[nmax+1]/(2m-1)!! is larger than 1/tiny.
 *             Tiny is defined below.
 *
 * If mode is equal to 2:
 *    Pl[l+1]*tiny/(2m-1)!!/(x^2-1)^(m/2)
 *    Ql[l+1]*(x^2-1)^(m/2)/(2m)!!/tiny
 *    nuevo    Maximum order of functions calculated when
 *             (nmax+1)Pl[nmax] is larger than 1/tiny.
 *             Tiny is defined below.
 */
void
prolate_harmonic(double x, int m, int nmax, int mode, std::vector<double> &Pl, std::vector<double> &Ql, int &nuevo)
{
    // dimension Pl(0:nmax+1),Ql(0:nmax+1)
    // required accuracy for the continued fraction (lentz-thompson)
    constexpr double eps = 1.0e-16;
    // Small parameter close to the underflow limit
    constexpr double tiny = 1.0e-290;
    constexpr auto huge = 1.0 / tiny;
    const auto sqrttiny = std::sqrt(tiny);
    if (x <= 1.0)
    {
        std::cout << "improper argument. x must be greater than 1";
        return;
    }
    const auto fl = m / 2.0;
    double ar;
    if (double(int(fl)) == fl)
        ar = 1.0;
    else
        ar = -1.0;
    const auto dzm = std::sqrt(x * x - 1.0);
    // We use the code if nmax is greater than or equal to 2
    auto nmaxp = nmax;
    if (nmax < 2)
        nmax = 2;
    // We store the values of Pl[0] and Pl[1]
    auto pl0 = 1.0;
    if (mode == 0)
    {
        if (m > 0)
        {
            int j = 1;
            for (int i = 1; i <= m; ++i)
            {
                pl0 *= dzm * (2 * m - j);
                j += 2;
            }
        }
        Pl[0] = pl0;
        Pl[1] = x * (2 * m + 1) * Pl[0];
    }
    else if (mode == 1)
    {
        if (std::abs(m * std::log(dzm)) > std::log(huge))
        {
            std::cout << "better try mode=2";
            return;
        }
        if (m > 0)
        {
            for (int i = 1; i <= m; ++i)
            {
                pl0 *= dzm;
            }
        }
        Pl[0] = pl0;
        Pl[1] = x * (2 * m + 1) * Pl[0];
    }
    else
    {
        Pl[0] = pl0 * tiny;
        Pl[1] = x * (2 * m + 1) * Pl[0];
    }

    // We use the recurrence relations
    int np = 1;
    while (np <= nmax && std::abs((np + 1.0) * Pl[np]) < huge)
    {
        Pl[np + 1] = ((2.0 * (np + m) + 1.0) * x * Pl[np] - (np + m + m) * Pl[np - 1]) / (np + 1.0);
        ++np;
    }
    nmax = np - 1;
    double factor;
    if (mode == 0)
    {
        factor = factco(nmax, Pl[nmax + 1], m, huge);
        while (factor == 0.0 && nmax > 0)
        {
            nmax = nmax - 3;
            factor = factco(nmax, Pl[nmax + 1], m, huge);
        }
        if (nmax <= 0)
        {
            std::cout << "try another m";
            return;
        }
    }
    else
    {
        factor = factconew(nmax, Pl[nmax + 1], m);
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
    a = -double(n + mm + 1) / double(n + 2 * m + mm);
    b = (1.0 + double(n + mm + 2) / double(2 * m + n + mm + 1)) * x;
    if (std::abs(delta - 1.0) > eps)
        goto _10;
    // We evaluate Ql[nmax+1] and Ql[nmax] using :
    // The Wronskian : W{Pl[nmax],Ql[nmax]} = 1 / (1 - x^2)
    // The known values of Pl[nmax+1] and Pl[nmax]
    // The value of h = Ql[nmax+1]/Ql[nmax]
    Ql[nmax] = ar * factor / (1.0 - fc * Pl[nmax] / Pl[nmax + 1]);
    Ql[nmax + 1] = Ql[nmax] * fc;
    for (int i = 1; i <= nmax; ++i)
    {
        np = nmax + 1 - i;
        n = np - 1;
        Ql[n] = ((2 * np + 2 * m + 1) * x * Ql[np] - (np + 1) * Ql[np + 1]) / double(np + m + m);
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
  std::vector<double> P(n_max + 2);
  std::vector<double> Q(n_max + 2);
  for (double x : {1.1, 1.5, 2.0, 5.0, 10.0})
  {
    prolate_harmonic(x, m, n_max, mode, P, Q, nuevo);
    std::cout << "\nx = " << x;
    std::cout << '\n';
    for (unsigned int n = 0; n <= n_max; ++n)
    {
      std::cout << ' ' << std::setw(12) << P[n];
    }
    std::cout << '\n';
    for (unsigned int n = 0; n <= n_max; ++n)
    {
      std::cout << ' ' << std::setw(12) << Q[n];
    }
  }
  std::cout << '\n';
}
