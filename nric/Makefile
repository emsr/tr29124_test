
SRCS = \
 bessel_general.c \
 bessel_integer.c \
 bessel_spherical.c \
 chebyshev.c \
 cholesky_decomp.c \
 dcomplex.c \
 diff.c \
 dlinear.c \
 dvector.c \
 fcomplex.c \
 fit.c \
 flinear.c \
 fourier.c \
 fvector.c \
 gamma.c \
 gauss_jordan.c \
 gauss_quad.c \
 int_cossin.c \
 int_dawson.c \
 int_elliptic.c \
 int_exp.c \
 int_fresnel.c \
 integration.c \
 interpolation.c \
 lu_decomp.c \
 memory.c \
 minimization.c \
 nric.c \
 ode.c \
 plot.c \
 poly.c \
 poly_chebyshev.c \
 poly_gegenbauer.c \
 poly_hermite.c \
 poly_jacobi.c \
 poly_laguerre.c \
 poly_legendre.c \
 qr_decomp.c \
 quadratic.c \
 roots.c \
 search_table.c \
 sorting.c \
 sph_harm.c \
 sv_decomp.c \
 tridiag.c util.c \
 wave_function.c

all: libnric.a graph integrate cheb matrix u vec lin

libnric.a: $(SRCS)
	gcc -std=c11 -static -o libnric.a  $(SRCS)

graph: graph.c plot.c gamma.c poly_*.c bessel_*.c util.c nric.c int_*.c dcomplex.c dvector.c fvector.c memory.c interpolation.c chebyshev.c poly.c integration.c gauss_quad.c
	gcc -std=c11 -o graph graph.c plot.c gamma.c poly_*.c bessel_*.c util.c nric.c int_*.c dcomplex.c dvector.c fvector.c memory.c interpolation.c chebyshev.c poly.c integration.c gauss_quad.c -lm

integrate: integrate.c integration.c chebyshev.c dvector.c fvector.c memory.c poly.c gauss_quad.c bessel_integer.c gamma.c interpolation.c plot.c nric.c
	gcc -std=c11 -o integrate integrate.c integration.c chebyshev.c dvector.c fvector.c memory.c poly.c gauss_quad.c bessel_integer.c gamma.c interpolation.c plot.c nric.c

cheb: cheb.c chebyshev.c dvector.c fvector.c memory.c poly.c nric.c 
	gcc -std=c11 -o cheb cheb.c chebyshev.c dvector.c fvector.c memory.c poly.c nric.c

matrix: matrix.c gauss_jordan.c lu_decomp.c qr_decomp.c sv_decomp.c cholesky_decomp.c dvector.c fvector.c memory.c nric.c dlinear.c flinear.c util.c
	gcc -std=c11 -o matrix matrix.c gauss_jordan.c lu_decomp.c qr_decomp.c sv_decomp.c cholesky_decomp.c dvector.c fvector.c memory.c nric.c dlinear.c flinear.c util.c

u: u.c util.c nric.c
	gcc -std=c11 -o u u.c util.c nric.c

vec: vec.c dvector.c fvector.c memory.c nric.c
	gcc -std=c11 -o vec vec.c dvector.c fvector.c memory.c nric.c

lin: lin.c dlinear.c dvector.c fvector.c memory.c nric.c
	gcc -std=c11 -o lin lin.c dlinear.c dvector.c fvector.c memory.c nric.c
