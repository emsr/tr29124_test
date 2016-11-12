
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

INC_DIR = .
OBJ_DIR = obj

OBJS = \
  $(OBJ_DIR)/drand.o \
  $(OBJ_DIR)/incbif.o \
  $(OBJ_DIR)/tandgf.o \
  $(OBJ_DIR)/ivf.o \
  $(OBJ_DIR)/j0f.o \
  $(OBJ_DIR)/expx2f.o \
  $(OBJ_DIR)/ellpjf.o \
  $(OBJ_DIR)/rgammaf.o \
  $(OBJ_DIR)/constf.o \
  $(OBJ_DIR)/gammaf.o \
  $(OBJ_DIR)/polevlf.o \
  $(OBJ_DIR)/powif.o \
  $(OBJ_DIR)/cmplxf.o \
  $(OBJ_DIR)/ellpef.o \
  $(OBJ_DIR)/expf.o \
  $(OBJ_DIR)/dawsnf.o \
  $(OBJ_DIR)/atanf.o \
  $(OBJ_DIR)/incbetf.o \
  $(OBJ_DIR)/igamf.o \
  $(OBJ_DIR)/sqrtf.o \
  $(OBJ_DIR)/ndtrf.o \
  $(OBJ_DIR)/psif.o \
  $(OBJ_DIR)/k0f.o \
  $(OBJ_DIR)/fresnlf.o \
  $(OBJ_DIR)/shichif.o \
  $(OBJ_DIR)/jvf.o \
  $(OBJ_DIR)/gdtrf.o \
  $(OBJ_DIR)/floorf.o \
  $(OBJ_DIR)/i0f.o \
  $(OBJ_DIR)/bdtrf.o \
  $(OBJ_DIR)/asinhf.o \
  $(OBJ_DIR)/logf.o \
  $(OBJ_DIR)/spencef.o \
  $(OBJ_DIR)/ynf.o \
  $(OBJ_DIR)/stdtrf.o \
  $(OBJ_DIR)/sindgf.o \
  $(OBJ_DIR)/hyp2f1f.o \
  $(OBJ_DIR)/cbrtf.o \
  $(OBJ_DIR)/powf.o \
  $(OBJ_DIR)/i1f.o \
  $(OBJ_DIR)/ellief.o \
  $(OBJ_DIR)/betaf.o \
  $(OBJ_DIR)/ndtrif.o \
  $(OBJ_DIR)/atanhf.o \
  $(OBJ_DIR)/expnf.o \
  $(OBJ_DIR)/tanf.o \
  $(OBJ_DIR)/k1f.o \
  $(OBJ_DIR)/tanhf.o \
  $(OBJ_DIR)/airyf.o \
  $(OBJ_DIR)/ellikf.o \
  $(OBJ_DIR)/asinf.o \
  $(OBJ_DIR)/ellpkf.o \
  $(OBJ_DIR)/exp2f.o \
  $(OBJ_DIR)/zetacf.o \
  $(OBJ_DIR)/j1f.o \
  $(OBJ_DIR)/hypergf.o \
  $(OBJ_DIR)/pdtrf.o \
  $(OBJ_DIR)/fdtrf.o \
  $(OBJ_DIR)/polynf.o \
  $(OBJ_DIR)/log2f.o \
  $(OBJ_DIR)/struvef.o \
  $(OBJ_DIR)/igamif.o \
  $(OBJ_DIR)/clogf.o \
  $(OBJ_DIR)/chdtrf.o \
  $(OBJ_DIR)/knf.o \
  $(OBJ_DIR)/facf.o \
  $(OBJ_DIR)/mtherr.o \
  $(OBJ_DIR)/zetaf.o \
  $(OBJ_DIR)/sicif.o \
  $(OBJ_DIR)/setprec.o \
  $(OBJ_DIR)/jnf.o \
  $(OBJ_DIR)/nbdtrf.o \
  $(OBJ_DIR)/sinf.o \
  $(OBJ_DIR)/exp10f.o \
  $(OBJ_DIR)/sinhf.o \
  $(OBJ_DIR)/acoshf.o \
  $(OBJ_DIR)/mtstf.o \
  $(OBJ_DIR)/log10f.o \
  $(OBJ_DIR)/chbevlf.o \
  $(OBJ_DIR)/coshf.o

libcephes_single.so: $(OBJS)
	$(CC) -shared -o libcephes_single.so $(OBJS) -L.. -lcephes_cmath
	cp libcephes_single.so libcephes_single.dll
	cp libcephes_single.so libcephes_single.dll ..

# /([a-z0-9]*)\.o
# /\1.o: \1.c\n\t\$\(CC\) -c -fPIC -I\$\(INC_DIR\) -o \$\(OBJ_DIR\)/\1.o \1.c\n

$(OBJ_DIR)/drand.o: drand.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/drand.o drand.c

$(OBJ_DIR)/incbif.o: incbif.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/incbif.o incbif.c

$(OBJ_DIR)/tandgf.o: tandgf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/tandgf.o tandgf.c

$(OBJ_DIR)/ivf.o: ivf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ivf.o ivf.c

$(OBJ_DIR)/j0f.o: j0f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/j0f.o j0f.c

$(OBJ_DIR)/expx2f.o: expx2f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/expx2f.o expx2f.c

$(OBJ_DIR)/ellpjf.o: ellpjf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellpjf.o ellpjf.c

$(OBJ_DIR)/rgammaf.o: rgammaf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/rgammaf.o rgammaf.c

$(OBJ_DIR)/constf.o: constf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/constf.o constf.c

$(OBJ_DIR)/gammaf.o: gammaf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/gammaf.o gammaf.c

$(OBJ_DIR)/polevlf.o: polevlf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polevlf.o polevlf.c

$(OBJ_DIR)/powif.o: powif.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/powif.o powif.c

$(OBJ_DIR)/cmplxf.o: cmplxf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cmplxf.o cmplxf.c

$(OBJ_DIR)/ellpef.o: ellpef.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellpef.o ellpef.c

$(OBJ_DIR)/expf.o: expf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/expf.o expf.c

$(OBJ_DIR)/dawsnf.o: dawsnf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/dawsnf.o dawsnf.c

$(OBJ_DIR)/atanf.o: atanf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/atanf.o atanf.c

$(OBJ_DIR)/incbetf.o: incbetf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/incbetf.o incbetf.c

$(OBJ_DIR)/igamf.o: igamf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/igamf.o igamf.c

$(OBJ_DIR)/sqrtf.o: sqrtf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sqrtf.o sqrtf.c

$(OBJ_DIR)/ndtrf.o: ndtrf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ndtrf.o ndtrf.c

$(OBJ_DIR)/psif.o: psif.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/psif.o psif.c

$(OBJ_DIR)/k0f.o: k0f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/k0f.o k0f.c

$(OBJ_DIR)/fresnlf.o: fresnlf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/fresnlf.o fresnlf.c

$(OBJ_DIR)/shichif.o: shichif.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/shichif.o shichif.c

$(OBJ_DIR)/jvf.o: jvf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/jvf.o jvf.c

$(OBJ_DIR)/gdtrf.o: gdtrf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/gdtrf.o gdtrf.c

$(OBJ_DIR)/floorf.o: floorf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/floorf.o floorf.c

$(OBJ_DIR)/i0f.o: i0f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/i0f.o i0f.c

$(OBJ_DIR)/bdtrf.o: bdtrf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/bdtrf.o bdtrf.c

$(OBJ_DIR)/asinhf.o: asinhf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/asinhf.o asinhf.c

$(OBJ_DIR)/logf.o: logf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/logf.o logf.c

$(OBJ_DIR)/spencef.o: spencef.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/spencef.o spencef.c

$(OBJ_DIR)/ynf.o: ynf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ynf.o ynf.c

$(OBJ_DIR)/stdtrf.o: stdtrf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/stdtrf.o stdtrf.c

$(OBJ_DIR)/sindgf.o: sindgf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sindgf.o sindgf.c

$(OBJ_DIR)/hyp2f1f.o: hyp2f1f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/hyp2f1f.o hyp2f1f.c

$(OBJ_DIR)/cbrtf.o: cbrtf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cbrtf.o cbrtf.c

$(OBJ_DIR)/powf.o: powf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/powf.o powf.c

$(OBJ_DIR)/i1f.o: i1f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/i1f.o i1f.c

$(OBJ_DIR)/ellief.o: ellief.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellief.o ellief.c

$(OBJ_DIR)/betaf.o: betaf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/betaf.o betaf.c

$(OBJ_DIR)/ndtrif.o: ndtrif.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ndtrif.o ndtrif.c

$(OBJ_DIR)/atanhf.o: atanhf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/atanhf.o atanhf.c

$(OBJ_DIR)/expnf.o: expnf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/expnf.o expnf.c

$(OBJ_DIR)/tanf.o: tanf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/tanf.o tanf.c

$(OBJ_DIR)/k1f.o: k1f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/k1f.o k1f.c

$(OBJ_DIR)/tanhf.o: tanhf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/tanhf.o tanhf.c

$(OBJ_DIR)/airyf.o: airyf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/airyf.o airyf.c

$(OBJ_DIR)/ellikf.o: ellikf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellikf.o ellikf.c

$(OBJ_DIR)/asinf.o: asinf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/asinf.o asinf.c

$(OBJ_DIR)/ellpkf.o: ellpkf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellpkf.o ellpkf.c

$(OBJ_DIR)/exp2f.o: exp2f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/exp2f.o exp2f.c

$(OBJ_DIR)/zetacf.o: zetacf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/zetacf.o zetacf.c

$(OBJ_DIR)/j1f.o: j1f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/j1f.o j1f.c

$(OBJ_DIR)/hypergf.o: hypergf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/hypergf.o hypergf.c

$(OBJ_DIR)/pdtrf.o: pdtrf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/pdtrf.o pdtrf.c

$(OBJ_DIR)/fdtrf.o: fdtrf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/fdtrf.o fdtrf.c

$(OBJ_DIR)/polynf.o: polynf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polynf.o polynf.c

$(OBJ_DIR)/log2f.o: log2f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/log2f.o log2f.c

$(OBJ_DIR)/struvef.o: struvef.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/struvef.o struvef.c

$(OBJ_DIR)/igamif.o: igamif.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/igamif.o igamif.c

$(OBJ_DIR)/clogf.o: clogf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/clogf.o clogf.c

$(OBJ_DIR)/chdtrf.o: chdtrf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/chdtrf.o chdtrf.c

$(OBJ_DIR)/knf.o: knf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/knf.o knf.c

$(OBJ_DIR)/facf.o: facf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/facf.o facf.c

$(OBJ_DIR)/mtherr.o: mtherr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mtherr.o mtherr.c

$(OBJ_DIR)/zetaf.o: zetaf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/zetaf.o zetaf.c

$(OBJ_DIR)/sicif.o: sicif.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sicif.o sicif.c

$(OBJ_DIR)/setprec.o: setprec.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/setprec.o setprec.c

$(OBJ_DIR)/jnf.o: jnf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/jnf.o jnf.c

$(OBJ_DIR)/nbdtrf.o: nbdtrf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/nbdtrf.o nbdtrf.c

$(OBJ_DIR)/sinf.o: sinf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sinf.o sinf.c

$(OBJ_DIR)/exp10f.o: exp10f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/exp10f.o exp10f.c

$(OBJ_DIR)/sinhf.o: sinhf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sinhf.o sinhf.c

$(OBJ_DIR)/acoshf.o: acoshf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/acoshf.o acoshf.c

$(OBJ_DIR)/mtstf.o: mtstf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mtstf.o mtstf.c

$(OBJ_DIR)/log10f.o: log10f.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/log10f.o log10f.c

$(OBJ_DIR)/chbevlf.o: chbevlf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/chbevlf.o chbevlf.c

$(OBJ_DIR)/coshf.o: coshf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/coshf.o coshf.c


clean:
	rm -f $(OBJS)


