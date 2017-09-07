
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

OBJ_DIR = obj

OBJS = \
 $(OBJ_DIR)/bessel_general.o \
 $(OBJ_DIR)/bessel_integer.o \
 $(OBJ_DIR)/bessel_spherical.o \
 $(OBJ_DIR)/chebyshev.o \
 $(OBJ_DIR)/cholesky_decomp.o \
 $(OBJ_DIR)/dcomplex.o \
 $(OBJ_DIR)/diff.o \
 $(OBJ_DIR)/dlinear.o \
 $(OBJ_DIR)/dvector.o \
 $(OBJ_DIR)/fcomplex.o \
 $(OBJ_DIR)/fit.o \
 $(OBJ_DIR)/flinear.o \
 $(OBJ_DIR)/fourier.o \
 $(OBJ_DIR)/fvector.o \
 $(OBJ_DIR)/gamma.o \
 $(OBJ_DIR)/gauss_jordan.o \
 $(OBJ_DIR)/gauss_quad.o \
 $(OBJ_DIR)/int_cossin.o \
 $(OBJ_DIR)/int_dawson.o \
 $(OBJ_DIR)/int_elliptic.o \
 $(OBJ_DIR)/int_exp.o \
 $(OBJ_DIR)/int_fresnel.o \
 $(OBJ_DIR)/integration.o \
 $(OBJ_DIR)/interpolation.o \
 $(OBJ_DIR)/lu_decomp.o \
 $(OBJ_DIR)/memory.o \
 $(OBJ_DIR)/minimization.o \
 $(OBJ_DIR)/nric.o \
 $(OBJ_DIR)/ode.o \
 $(OBJ_DIR)/plot.o \
 $(OBJ_DIR)/poly.o \
 $(OBJ_DIR)/poly_chebyshev.o \
 $(OBJ_DIR)/poly_gegenbauer.o \
 $(OBJ_DIR)/poly_hermite.o \
 $(OBJ_DIR)/poly_jacobi.o \
 $(OBJ_DIR)/poly_laguerre.o \
 $(OBJ_DIR)/poly_legendre.o \
 $(OBJ_DIR)/qr_decomp.o \
 $(OBJ_DIR)/quadratic.o \
 $(OBJ_DIR)/roots.o \
 $(OBJ_DIR)/search_table.o \
 $(OBJ_DIR)/sorting.o \
 $(OBJ_DIR)/sph_harm.o \
 $(OBJ_DIR)/sv_decomp.o \
 $(OBJ_DIR)/tridiag.o util.o \
 $(OBJ_DIR)/wave_function.o

all: libnric.a graph integrate cheb matrix u vec lin

libnric.a: $(OBJ_DIR) $(OBJS)
	ar cr libnric.a  $(OBJS)

graph: graph.c plot.c gamma.c poly_*.c bessel_*.c util.c nric.c int_*.c dcomplex.c dvector.c fvector.c memory.c interpolation.c chebyshev.c poly.c integration.c gauss_quad.c
	gcc -std=c11 -o graph graph.c plot.c gamma.c poly_*.c bessel_*.c util.c nric.c int_*.c dcomplex.c dvector.c fvector.c memory.c interpolation.c chebyshev.c poly.c integration.c gauss_quad.c -lm

integrate: integrate.c integration.c chebyshev.c dvector.c fvector.c memory.c poly.c gauss_quad.c bessel_integer.c gamma.c interpolation.c plot.c nric.c
	g++ -std=c++11 -o integrate integrate.c integration.c chebyshev.c dvector.c fvector.c memory.c poly.c gauss_quad.c bessel_integer.c gamma.c interpolation.c plot.c nric.c -lm

cheb: cheb.c chebyshev.c dvector.c fvector.c memory.c poly.c nric.c 
	gcc -std=c11 -o cheb cheb.c chebyshev.c dvector.c fvector.c memory.c poly.c nric.c -lm

matrix: matrix.c gauss_jordan.c lu_decomp.c qr_decomp.c sv_decomp.c cholesky_decomp.c dvector.c fvector.c memory.c nric.c dlinear.c flinear.c util.c
	gcc -std=c11 -o matrix matrix.c gauss_jordan.c lu_decomp.c qr_decomp.c sv_decomp.c cholesky_decomp.c dvector.c fvector.c memory.c nric.c dlinear.c flinear.c util.c -lm

u: u.c util.c nric.c
	gcc -std=c11 -o u u.c util.c nric.c -lm

vec: vec.c dvector.c fvector.c memory.c nric.c
	gcc -std=c11 -o vec vec.c dvector.c fvector.c memory.c nric.c -lm

lin: lin.c dlinear.c dvector.c fvector.c memory.c nric.c
	gcc -std=c11 -o lin lin.c dlinear.c dvector.c fvector.c memory.c nric.c -lm


$(OBJ_DIR)/bessel_general.o: bessel_general.c
	gcc -std=c11 -c -o $(OBJ_DIR)/bessel_general.o bessel_general.c

$(OBJ_DIR)/bessel_integer.o: bessel_integer.c
	gcc -std=c11 -c -o $(OBJ_DIR)/bessel_integer.o bessel_integer.c

$(OBJ_DIR)/bessel_spherical.o: bessel_spherical.c
	gcc -std=c11 -c -o $(OBJ_DIR)/bessel_spherical.o bessel_spherical.c

$(OBJ_DIR)/chebyshev.o: chebyshev.c
	gcc -std=c11 -c -o $(OBJ_DIR)/chebyshev.o chebyshev.c

$(OBJ_DIR)/cholesky_decomp.o: cholesky_decomp.c
	gcc -std=c11 -c -o $(OBJ_DIR)/cholesky_decomp.o cholesky_decomp.c

$(OBJ_DIR)/dcomplex.o: dcomplex.c
	gcc -std=c11 -c -o $(OBJ_DIR)/dcomplex.o dcomplex.c

$(OBJ_DIR)/diff.o: diff.c
	gcc -std=c11 -c -o $(OBJ_DIR)/diff.o diff.c

$(OBJ_DIR)/dlinear.o: dlinear.c
	gcc -std=c11 -c -o $(OBJ_DIR)/dlinear.o dlinear.c

$(OBJ_DIR)/dvector.o: dvector.c
	gcc -std=c11 -c -o $(OBJ_DIR)/dvector.o dvector.c

$(OBJ_DIR)/fcomplex.o: fcomplex.c
	gcc -std=c11 -c -o $(OBJ_DIR)/fcomplex.o fcomplex.c

$(OBJ_DIR)/fit.o: fit.c
	gcc -std=c11 -c -o $(OBJ_DIR)/fit.o fit.c

$(OBJ_DIR)/flinear.o: flinear.c
	gcc -std=c11 -c -o $(OBJ_DIR)/flinear.o flinear.c

$(OBJ_DIR)/fourier.o: fourier.c
	gcc -std=c11 -c -o $(OBJ_DIR)/fourier.o fourier.c

$(OBJ_DIR)/fvector.o: fvector.c
	gcc -std=c11 -c -o $(OBJ_DIR)/fvector.o fvector.c

$(OBJ_DIR)/gamma.o: gamma.c
	gcc -std=c11 -c -o $(OBJ_DIR)/gamma.o gamma.c

$(OBJ_DIR)/gauss_jordan.o: gauss_jordan.c
	gcc -std=c11 -c -o $(OBJ_DIR)/gauss_jordan.o gauss_jordan.c

$(OBJ_DIR)/gauss_quad.o: gauss_quad.c
	gcc -std=c11 -c -o $(OBJ_DIR)/gauss_quad.o gauss_quad.c

$(OBJ_DIR)/int_cossin.o: int_cossin.c
	gcc -std=c11 -c -o $(OBJ_DIR)/int_cossin.o int_cossin.c

$(OBJ_DIR)/int_dawson.o: int_dawson.c
	gcc -std=c11 -c -o $(OBJ_DIR)/int_dawson.o int_dawson.c

$(OBJ_DIR)/int_elliptic.o: int_elliptic.c
	gcc -std=c11 -c -o $(OBJ_DIR)/int_elliptic.o int_elliptic.c

$(OBJ_DIR)/int_exp.o: int_exp.c
	gcc -std=c11 -c -o $(OBJ_DIR)/int_exp.o int_exp.c

$(OBJ_DIR)/int_fresnel.o: int_fresnel.c
	gcc -std=c11 -c -o $(OBJ_DIR)/int_fresnel.o int_fresnel.c

$(OBJ_DIR)/integration.o: integration.c
	gcc -std=c11 -c -o $(OBJ_DIR)/integration.o integration.c

$(OBJ_DIR)/interpolation.o: interpolation.c
	gcc -std=c11 -c -o $(OBJ_DIR)/interpolation.o interpolation.c

$(OBJ_DIR)/lu_decomp.o: lu_decomp.c
	gcc -std=c11 -c -o $(OBJ_DIR)/lu_decomp.o lu_decomp.c

$(OBJ_DIR)/memory.o: memory.c
	gcc -std=c11 -c -o $(OBJ_DIR)/memory.o memory.c

$(OBJ_DIR)/minimization.o: minimization.c
	gcc -std=c11 -c -o $(OBJ_DIR)/minimization.o minimization.c

$(OBJ_DIR)/nric.o: nric.c
	gcc -std=c11 -c -o $(OBJ_DIR)/nric.o nric.c

$(OBJ_DIR)/ode.o: ode.c
	gcc -std=c11 -c -o $(OBJ_DIR)/ode.o ode.c

$(OBJ_DIR)/plot.o: plot.c
	gcc -std=c11 -c -o $(OBJ_DIR)/plot.o plot.c

$(OBJ_DIR)/poly.o: poly.c
	gcc -std=c11 -c -o $(OBJ_DIR)/poly.o poly.c

$(OBJ_DIR)/poly_chebyshev.o: poly_chebyshev.c
	gcc -std=c11 -c -o $(OBJ_DIR)/poly_chebyshev.o poly_chebyshev.c

$(OBJ_DIR)/poly_gegenbauer.o: poly_gegenbauer.c
	gcc -std=c11 -c -o $(OBJ_DIR)/poly_gegenbauer.o poly_gegenbauer.c

$(OBJ_DIR)/poly_hermite.o: poly_hermite.c
	gcc -std=c11 -c -o $(OBJ_DIR)/poly_hermite.o poly_hermite.c

$(OBJ_DIR)/poly_jacobi.o: poly_jacobi.c
	gcc -std=c11 -c -o $(OBJ_DIR)/poly_jacobi.o poly_jacobi.c

$(OBJ_DIR)/poly_laguerre.o: poly_laguerre.c
	gcc -std=c11 -c -o $(OBJ_DIR)/poly_laguerre.o poly_laguerre.c

$(OBJ_DIR)/poly_legendre.o: poly_legendre.c
	gcc -std=c11 -c -o $(OBJ_DIR)/poly_legendre.o poly_legendre.c

$(OBJ_DIR)/qr_decomp.o: qr_decomp.c
	gcc -std=c11 -c -o $(OBJ_DIR)/qr_decomp.o qr_decomp.c

$(OBJ_DIR)/quadratic.o: quadratic.c
	gcc -std=c11 -c -o $(OBJ_DIR)/quadratic.o quadratic.c

$(OBJ_DIR)/roots.o: roots.c
	gcc -std=c11 -c -o $(OBJ_DIR)/roots.o roots.c

$(OBJ_DIR)/search_table.o: search_table.c
	gcc -std=c11 -c -o $(OBJ_DIR)/search_table.o search_table.c

$(OBJ_DIR)/sorting.o: sorting.c
	gcc -std=c11 -c -o $(OBJ_DIR)/sorting.o sorting.c

$(OBJ_DIR)/sph_harm.o: sph_harm.c
	gcc -std=c11 -c -o $(OBJ_DIR)/sph_harm.o sph_harm.c

$(OBJ_DIR)/sv_decomp.o: sv_decomp.c
	gcc -std=c11 -c -o $(OBJ_DIR)/sv_decomp.o sv_decomp.c

$(OBJ_DIR)/tridiag.o: tridiag.c
	gcc -std=c11 -c -o $(OBJ_DIR)/tridiag.o tridiag.c

$(OBJ_DIR)/util.o: util.c
	gcc -std=c11 -c -o $(OBJ_DIR)/util.o util.c

$(OBJ_DIR)/wave_function.o: wave_function.c
	gcc -std=c11 -c -o $(OBJ_DIR)/wave_function.o wave_function.c



$(OBJ_DIR): $(OUT_DIR)
	if test ! -d $(OBJ_DIR); then \
	  mkdir $(OBJ_DIR); \
	fi
