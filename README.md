This repository contains work toward Mathematical Special Functions in C++.
This is targeted towards (and licensed for) libstdc++ but I hope that
it could be generally useful and amusing.

This was originally in TR29124.  Hence the awkward name.

Now these functions have been accepted into C++17.

In addition to the special functions in C++17, this library adds several extensions:
* Hypergeometric functions,
* Carlson elliptic functions,
* Polylogarithm functions,
* Hankel functions,
* Statistical functions
* Quadrature rules,
,...

I strive for type genericity.  I want C++ numerics to follow
the containers + algorithms by having numeriic algorithms that will
work for any type for which numeric_limits, and the basic math functions
are available.  These functions have been tested with float, double, long double,
and __float128.  Some have been tested with mpreal and efforts are underway to
allow full multiprecision usage.

I generally strive towards accuracy first and speed second.
