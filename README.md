# TR29124 Mathematical Special Functions in C++

This repository contains work toward
[IS 29124 - Extensions to the C++ Library to Support Mathematical Special Functions]
(http://www.open-std.org/JTC1/SC22/WG21/docs/papers/2010/n3060.pdf)
hence the awkward name.

This library also contains work towards follow-up published proposals
for new special functions:
[A proposal to add special mathematical functions according to
the ISO/IEC 80000-2:2009 standard, Vincent Reverdy]
(http://open-std.org/JTC1/SC22/WG21/docs/papers/2013/n3494.pdf)

[A Proposal to add Mathematical Functions for Statistics
to the C++ Standard Library, Paul A Bristow]
(http://open-std.org/JTC1/SC22/WG21/docs/papers/2004/n1668.pdf)

[A proposal to add sincos to the standard library, Paul Dreik]
(http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/p0081r0.pdf)

This work is targeted towards (and licensed for) libstdc++ but I hope that
it is generally useful and amusing.

The functions in IS 29124 have been accepted into C++17.
See Section 26.9.5 Mathematical special functions [sf.cmath] in a recent draft.

In addition to the special functions in C++17, this library adds several extensions:
  * Hypergeometric functions
  * Carlson elliptic functions
  * Jacobi elliptic functions, amplitude, and nome
  * Polylogarithm functions
  * Airy functions
  * Hankel functions
  * Statistical distribution functions

Some capabilities have bee migrated into separate subprojects:
* Quadrature rules
* Polynomial classes and solvers

I strive for type genericity.  I want C++ numerics to follow
the containers + algorithms by having numeric algorithms that will
work for any type for which numeric_limits, and the basic math functions
are available.  These functions have been tested with float, double, long double,
and \_\_float128.  Some have been tested with mpreal and efforts are underway to
allow full multiprecision usage.

I generally strive towards accuracy first and speed second.

## Documentation

Check out the [API documentation](https://emsr.github.io/tr29124_test/html/index.html).

## Building

Because there are several submodules, when you first clone this package you'll neeed to go
into the main project directory and update the submodules:
  * git submodule update --init --recursive
From then on:
  * git submodule update --recursive --remote

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/bee73172954646a0bf287bdcdd269553)](https://app.codacy.com/app/emsr/tr29124_test?utm_source=github.com&utm_medium=referral&utm_content=emsr/tr29124_test&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.org/emsr/tr29124_test.svg?branch=master)](https://travis-ci.org/emsr/tr29124_test)
[![Coverage Status](https://codecov.io/gh/emsr/tr29124_test/branch/master/graph/badge.svg)](https://codecov.io/gh/emsr/tr29124_test)
