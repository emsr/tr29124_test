//  cgamma.cpp -- Complex gamma function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns gamma function or log(gamma) for complex argument 'z'.
//
//  OPT value       function
//  ---------       --------
//      0           complex gamma
//      1           complex log(gamma)
//
//  Returns (1e308,0) if the real part of the argument is a negative integer
//  or 0 or exceeds 171.
//

#include <complex>

template<typename Tp>
  complex<Tp>
  cgamma(const std::complex<Tp>& z, int OPT)
  {
    complex<Tp> z0, z1;

    static Tp a[] = {
        8.333333333333333e-02,
       -2.777777777777778e-03,
        7.936507936507937e-04,
       -5.952380952380952e-04,
        8.417508417508418e-04,
       -1.917526917526918e-03,
        6.410256410256410e-03,
       -2.955065359477124e-02,
        1.796443723688307e-01,
       -1.39243221690590
    };

    Tp x = std::real(z);
    Tp y = std::imag(z);
    if (x > 171)
      return std::complex<Tp>(1e308, Tp(0));
    if ((y == Tp(0)) && (x == (int)x) && (x <= Tp(0)))
      return std::complex<Tp>(1e308, Tp(0));
    else if (x < Tp(0))
      {
        x1 = x;
        y1 = y;
        x = -x;
        y = -y;
      }
    int na = 0;
    Tp x0 = x;
    if (x <= 7.0) 
      {
        na = (int)(7.0 - x);
        x0 = x + na;
      }
    Tp q1 = std::sqrt(x0 * x0 + y * y);
    Tp th = std::atan2(y, x0);
    Tp gr = (x0 - 0.5) * std::log(q1) - th * y - x0 + 0.5 * std::log(2.0 * M_PI);
    Tp gi = th * (x0 - 0.5) + y * std::log(q1) - y;
    for (std::size_t k = 0; k < 10; ++k)
      {
        Tp t = std::pow(q1, -Tp(2 * k + 1));
        gr += a[k] * t * std::cos(Tp(2 * k + 1) * th);
        gi -= a[k] * t * std::sin(Tp(2 * k + 1) * th);
      }
    if (x <= 7.0)
      {
        Tp gr1 = Tp(0);
        Tp gi1 = Tp(0);
        for (int j = 0; j < na; ++j)
          {
            gr1 += (0.5 * std::log((x + j) * (x + j) + y * y));
            gi1 += std::atan2(y, x + j);
          }
        gr -= gr1;
        gi -= gi1;
      }
    if (x1 <= Tp(0))
      {
        q1 = std::sqrt(x * x + y * y);
        Tp th1 = std::atan2(y, x);
        Tp sr = -std::sin(M_PI * x) * std::cosh(M_PI * y);
        Tp si = -std::cos(M_PI * x) * std::sinh(M_PI * y);
        Tp q2 = std::sqrt(sr * sr + si * si);
        Tp th2 = std::atan2(si, sr);
        if (sr < Tp(0))
          th2 += M_PI;
        gr = std::log(M_PI / (q1 * q2)) - gr;
        gi = -th1 - th2 - gi;
        x = x1;
        y = y1;
      }
    if (OPT == 0)
      {
        Tp g0 = std::exp(gr);
        gr = g0 * std::cos(gi);
        gi = g0 * std::sin(gi);
      }

    return std::complex<Tp>(gr, gi);
  }
