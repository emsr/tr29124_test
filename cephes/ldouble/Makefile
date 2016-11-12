
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

INC_DIR = .
OBJ_DIR = obj

SRCS = \
  acoshl.c \
  asinl.c \
  asinhl.c \
  atanl.c \
  atanhl.c \
  bdtrl.c \
  btdtrl.c \
  cbrtl.c \
  clogl.c \
  cmplxl.c \
  coshl.c \
  econst.c \
  elliel.c \
  ellikl.c \
  ellpel.c \
  ellpjl.c \
  ellpkl.c \
  expl.c \
  exp10l.c \
  exp2l.c \
  floorl.c \
  gammal.c \
  gdtrl.c \
  gelsl.c \
  ieee.c \
  igamil.c \
  igaml.c \
  incbetl.c \
  incbil.c \
  isnanl.c \
  j0l.c \
  j1l.c \
  jnl.c \
  logl.c \
  log10l.c \
  log2l.c \
  mod2pil.c \
  mtstl.c \
  nbdtrl.c \
  ndtril.c \
  ndtrl.c \
  pdtrl.c \
  polevll.c \
  powl.c \
  powil.c \
  sinl.c \
  sinhl.c \
  sqrtl.c \
  stdtrl.c \
  tanl.c \
  tanhl.c \
  unityl.c \
  ynl.c

OBJS = \
  $(OBJ_DIR)/acoshl.o \
  $(OBJ_DIR)/asinl.o \
  $(OBJ_DIR)/asinhl.o \
  $(OBJ_DIR)/atanl.o \
  $(OBJ_DIR)/atanhl.o \
  $(OBJ_DIR)/bdtrl.o \
  $(OBJ_DIR)/btdtrl.o \
  $(OBJ_DIR)/cbrtl.o \
  $(OBJ_DIR)/clogl.o \
  $(OBJ_DIR)/cmplxl.o \
  $(OBJ_DIR)/coshl.o \
  $(OBJ_DIR)/econst.o \
  $(OBJ_DIR)/elliel.o \
  $(OBJ_DIR)/ellikl.o \
  $(OBJ_DIR)/ellpel.o \
  $(OBJ_DIR)/ellpjl.o \
  $(OBJ_DIR)/ellpkl.o \
  $(OBJ_DIR)/expl.o \
  $(OBJ_DIR)/expx2l.o \
  $(OBJ_DIR)/exp10l.o \
  $(OBJ_DIR)/exp2l.o \
  $(OBJ_DIR)/floorl.o \
  $(OBJ_DIR)/gammal.o \
  $(OBJ_DIR)/gdtrl.o \
  $(OBJ_DIR)/gelsl.o \
  $(OBJ_DIR)/ieee.o \
  $(OBJ_DIR)/igamil.o \
  $(OBJ_DIR)/igaml.o \
  $(OBJ_DIR)/incbetl.o \
  $(OBJ_DIR)/incbil.o \
  $(OBJ_DIR)/isnanl.o \
  $(OBJ_DIR)/j0l.o \
  $(OBJ_DIR)/j1l.o \
  $(OBJ_DIR)/jnl.o \
  $(OBJ_DIR)/logl.o \
  $(OBJ_DIR)/log10l.o \
  $(OBJ_DIR)/log2l.o \
  $(OBJ_DIR)/nbdtrl.o \
  $(OBJ_DIR)/ndtril.o \
  $(OBJ_DIR)/ndtrl.o \
  $(OBJ_DIR)/pdtrl.o \
  $(OBJ_DIR)/polevll.o \
  $(OBJ_DIR)/powl.o \
  $(OBJ_DIR)/powil.o \
  $(OBJ_DIR)/sinl.o \
  $(OBJ_DIR)/sinhl.o \
  $(OBJ_DIR)/sqrtl.o \
  $(OBJ_DIR)/stdtrl.o \
  $(OBJ_DIR)/tanl.o \
  $(OBJ_DIR)/tanhl.o \
  $(OBJ_DIR)/unityl.o \
  $(OBJ_DIR)/ynl.o

libcephes_ldouble.so: $(OBJS)
	$(CC) -shared -o libcephes_ldouble.so $(OBJS) -L.. -lcephes_cmath
	cp libcephes_ldouble.so libcephes_ldouble.dll
	cp libcephes_ldouble.so libcephes_ldouble.dll ..

# /([a-z0-9]*)\.o
# /\1.o: \1.c\n\t\$\(CC\) -c -fPIC -I\$\(INC_DIR\) -o \$\(OBJ_DIR\)/\1.o \1.c\n

$(OBJ_DIR)/acoshl.o: acoshl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/acoshl.o acoshl.c

$(OBJ_DIR)/asinl.o: asinl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/asinl.o asinl.c

$(OBJ_DIR)/asinhl.o: asinhl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/asinhl.o asinhl.c

$(OBJ_DIR)/atanl.o: atanl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/atanl.o atanl.c

$(OBJ_DIR)/atanhl.o: atanhl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/atanhl.o atanhl.c

$(OBJ_DIR)/bdtrl.o: bdtrl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/bdtrl.o bdtrl.c

$(OBJ_DIR)/btdtrl.o: btdtrl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/btdtrl.o btdtrl.c

$(OBJ_DIR)/cbrtl.o: cbrtl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cbrtl.o cbrtl.c

$(OBJ_DIR)/clogl.o: clogl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/clogl.o clogl.c

$(OBJ_DIR)/cmplxl.o: cmplxl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cmplxl.o cmplxl.c

$(OBJ_DIR)/coshl.o: coshl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/coshl.o coshl.c

$(OBJ_DIR)/econst.o: econst.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/econst.o econst.c

$(OBJ_DIR)/elliel.o: elliel.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/elliel.o elliel.c

$(OBJ_DIR)/ellikl.o: ellikl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellikl.o ellikl.c

$(OBJ_DIR)/ellpel.o: ellpel.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellpel.o ellpel.c

$(OBJ_DIR)/ellpjl.o: ellpjl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellpjl.o ellpjl.c

$(OBJ_DIR)/ellpkl.o: ellpkl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellpkl.o ellpkl.c

$(OBJ_DIR)/dtestvecl.o: dtestvecl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/dtestvecl.o dtestvecl.c

$(OBJ_DIR)/expl.o: expl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/expl.o expl.c

$(OBJ_DIR)/expx2l.o: expx2l.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/expx2l.o expx2l.c

$(OBJ_DIR)/exp10l.o: exp10l.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/exp10l.o exp10l.c

$(OBJ_DIR)/exp2l.o: exp2l.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/exp2l.o exp2l.c

$(OBJ_DIR)/floorl.o: floorl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/floorl.o floorl.c

$(OBJ_DIR)/gammal.o: gammal.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/gammal.o gammal.c

$(OBJ_DIR)/gdtrl.o: gdtrl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/gdtrl.o gdtrl.c

$(OBJ_DIR)/gelsl.o: gelsl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/gelsl.o gelsl.c

$(OBJ_DIR)/ieee.o: ieee.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ieee.o ieee.c

$(OBJ_DIR)/igamil.o: igamil.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/igamil.o igamil.c

$(OBJ_DIR)/igaml.o: igaml.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/igaml.o igaml.c

$(OBJ_DIR)/incbetl.o: incbetl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/incbetl.o incbetl.c

$(OBJ_DIR)/incbil.o: incbil.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/incbil.o incbil.c

$(OBJ_DIR)/isnanl.o: isnanl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/isnanl.o isnanl.c

$(OBJ_DIR)/j0l.o: j0l.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/j0l.o j0l.c

$(OBJ_DIR)/j1l.o: j1l.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/j1l.o j1l.c

$(OBJ_DIR)/jnl.o: jnl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/jnl.o jnl.c

$(OBJ_DIR)/logl.o: logl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/logl.o logl.c

$(OBJ_DIR)/log10l.o: log10l.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/log10l.o log10l.c

$(OBJ_DIR)/log2l.o: log2l.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/log2l.o log2l.c

$(OBJ_DIR)/mod2pil.o: mod2pil.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mod2pil.o mod2pil.c

$(OBJ_DIR)/mtstl.o: mtstl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mtstl.o mtstl.c

$(OBJ_DIR)/nbdtrl.o: nbdtrl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/nbdtrl.o nbdtrl.c

$(OBJ_DIR)/ndtril.o: ndtril.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ndtril.o ndtril.c

$(OBJ_DIR)/ndtrl.o: ndtrl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ndtrl.o ndtrl.c

$(OBJ_DIR)/pdtrl.o: pdtrl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/pdtrl.o pdtrl.c

$(OBJ_DIR)/polevll.o: polevll.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polevll.o polevll.c

$(OBJ_DIR)/powl.o: powl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/powl.o powl.c

$(OBJ_DIR)/powil.o: powil.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/powil.o powil.c

$(OBJ_DIR)/sinl.o: sinl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sinl.o sinl.c

$(OBJ_DIR)/sinhl.o: sinhl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sinhl.o sinhl.c

$(OBJ_DIR)/sqrtl.o: sqrtl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sqrtl.o sqrtl.c

$(OBJ_DIR)/stdtrl.o: stdtrl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/stdtrl.o stdtrl.c

$(OBJ_DIR)/tanl.o: tanl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/tanl.o tanl.c

$(OBJ_DIR)/tanhl.o: tanhl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/tanhl.o tanhl.c

$(OBJ_DIR)/unityl.o: unityl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/unityl.o unityl.c

$(OBJ_DIR)/ynl.o: ynl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ynl.o ynl.c


clean:
	rm -f $(OBJS)

