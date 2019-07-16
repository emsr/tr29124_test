Add a Doxy mainpage and a README.md

Why does the test print out values that look like long double rather than __float128 at te end?
It's the definition of pi_v as cast from a long double!  We need some kind of large thing that
can feed whatever bits people need.
