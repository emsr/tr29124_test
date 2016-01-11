all: graph integrate cheb matrix u vec lin

graph: graph.c plot.c gamma.c poly_*.c bessel_*.c util.c nric.c int_*.c dcomplex.c dvector.c fvector.c memory.c interpolation.c chebyshev.c poly.c integration.c gauss_quad.c
	gcc -o graph graph.c plot.c gamma.c poly_*.c bessel_*.c util.c nric.c int_*.c dcomplex.c dvector.c fvector.c memory.c interpolation.c chebyshev.c poly.c integration.c gauss_quad.c -lm

integrate: integrate.c integration.c chebyshev.c dvector.c fvector.c memory.c poly.c gauss_quad.c bessel_integer.c gamma.c interpolation.c plot.c nric.c
	gcc -o integrate integrate.c integration.c chebyshev.c dvector.c fvector.c memory.c poly.c gauss_quad.c bessel_integer.c gamma.c interpolation.c plot.c nric.c

cheb: cheb.c chebyshev.c dvector.c fvector.c memory.c poly.c nric.c 
	gcc -o cheb cheb.c chebyshev.c dvector.c fvector.c memory.c poly.c nric.c

matrix: matrix.c gauss_jordan.c lu_decomp.c qr_decomp.c sv_decomp.c cholesky_decomp.c dvector.c fvector.c memory.c nric.c dlinear.c flinear.c util.c
	gcc -o matrix matrix.c gauss_jordan.c lu_decomp.c qr_decomp.c sv_decomp.c cholesky_decomp.c dvector.c fvector.c memory.c nric.c dlinear.c flinear.c util.c

u: u.c util.c nric.c
	gcc -o u u.c util.c nric.c

vec: vec.c dvector.c fvector.c memory.c nric.c
	gcc -o vec vec.c dvector.c fvector.c memory.c nric.c

lin: lin.c dlinear.c dvector.c fvector.c memory.c nric.c
	gcc -o lin lin.c dlinear.c dvector.c fvector.c memory.c nric.c
