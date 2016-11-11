
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

INC_DIR = .
OBJ_DIR = obj

SRCS = \
  acosh.c \
  asin.c \
  asinh.c \
  atan.c \
  atanh.c \
  cbrt.c \
  chbevl.c \
  clog.c \
  cmplx.c \
  const.c \
  cosh.c \
  drand.c \
  dtestvec.c \
  exp.c \
  exp10.c \
  exp2.c \
  fabs.c \
  floor.c \
  isnan.c \
  log.c \
  log10.c \
  log2.c \
  mod2pi.c \
  mtherr.c \
  mtst.c \
  polevl.c \
  pow.c \
  powi.c \
  round.c \
  setprec.c \
  sin.c \
  sincos.c \
  sindg.c \
  sinh.c \
  sqrt.c \
  tan.c \
  tandg.c \
  tanh.c \
  unity.c

#  $(OBJ_DIR)/mtst.o
#  $(OBJ_DIR)/mod2pi.o
#  $(OBJ_DIR)/dtestvec.o
OBJS = \
  $(OBJ_DIR)/acosh.o \
  $(OBJ_DIR)/asin.o \
  $(OBJ_DIR)/asinh.o \
  $(OBJ_DIR)/atan.o \
  $(OBJ_DIR)/atanh.o \
  $(OBJ_DIR)/cbrt.o \
  $(OBJ_DIR)/chbevl.o \
  $(OBJ_DIR)/clog.o \
  $(OBJ_DIR)/cmplx.o \
  $(OBJ_DIR)/const.o \
  $(OBJ_DIR)/cosh.o \
  $(OBJ_DIR)/drand.o \
  $(OBJ_DIR)/exp.o \
  $(OBJ_DIR)/exp10.o \
  $(OBJ_DIR)/exp2.o \
  $(OBJ_DIR)/fabs.o \
  $(OBJ_DIR)/floor.o \
  $(OBJ_DIR)/isnan.o \
  $(OBJ_DIR)/log.o \
  $(OBJ_DIR)/log10.o \
  $(OBJ_DIR)/log2.o \
  $(OBJ_DIR)/mtherr.o \
  $(OBJ_DIR)/polevl.o \
  $(OBJ_DIR)/pow.o \
  $(OBJ_DIR)/powi.o \
  $(OBJ_DIR)/round.o \
  $(OBJ_DIR)/setprec.o \
  $(OBJ_DIR)/sin.o \
  $(OBJ_DIR)/sincos.o \
  $(OBJ_DIR)/sindg.o \
  $(OBJ_DIR)/sinh.o \
  $(OBJ_DIR)/sqrt.o \
  $(OBJ_DIR)/tan.o \
  $(OBJ_DIR)/tandg.o \
  $(OBJ_DIR)/tanh.o \
  $(OBJ_DIR)/unity.o

libcephes_cmath.so: $(OBJS)
	$(CC) -shared -o libcephes_cmath.so $(OBJS)
	cp libcephes_cmath.so libcephes_cmath.dll
	cp libcephes_cmath.so libcephes_cmath.dll ..

$(OBJ_DIR)/acosh.o: acosh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/acosh.o acosh.c

$(OBJ_DIR)/asin.o: asin.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/asin.o asin.c

$(OBJ_DIR)/asinh.o: asinh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/asinh.o asinh.c

$(OBJ_DIR)/atan.o: atan.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/atan.o atan.c

$(OBJ_DIR)/atanh.o: atanh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/atanh.o atanh.c

$(OBJ_DIR)/cbrt.o: cbrt.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cbrt.o cbrt.c

$(OBJ_DIR)/chbevl.o: chbevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/chbevl.o chbevl.c

$(OBJ_DIR)/clog.o: clog.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/clog.o clog.c

$(OBJ_DIR)/cmplx.o: cmplx.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cmplx.o cmplx.c

$(OBJ_DIR)/const.o: const.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/const.o const.c

$(OBJ_DIR)/cosh.o: cosh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cosh.o cosh.c

$(OBJ_DIR)/drand.o: drand.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/drand.o drand.c

$(OBJ_DIR)/dtestvec.o: dtestvec.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/dtestvec.o dtestvec.c

$(OBJ_DIR)/exp.o: exp.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/exp.o exp.c

$(OBJ_DIR)/exp10.o: exp10.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/exp10.o exp10.c

$(OBJ_DIR)/exp2.o: exp2.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/exp2.o exp2.c

$(OBJ_DIR)/fabs.o: fabs.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/fabs.o fabs.c

$(OBJ_DIR)/floor.o: floor.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/floor.o floor.c

$(OBJ_DIR)/isnan.o: isnan.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/isnan.o isnan.c

$(OBJ_DIR)/log.o: log.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/log.o log.c

$(OBJ_DIR)/log10.o: log10.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/log10.o log10.c

$(OBJ_DIR)/log2.o: log2.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/log2.o log2.c

$(OBJ_DIR)/mod2pi.o: mod2pi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mod2pi.o mod2pi.c

$(OBJ_DIR)/mtherr.o: mtherr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mtherr.o mtherr.c

$(OBJ_DIR)/mtst.o: mtst.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mtst.o mtst.c

$(OBJ_DIR)/polevl.o: polevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polevl.o polevl.c

$(OBJ_DIR)/pow.o: pow.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/pow.o pow.c

$(OBJ_DIR)/powi.o: powi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/powi.o powi.c

$(OBJ_DIR)/round.o: round.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/round.o round.c

$(OBJ_DIR)/setprec.o: setprec.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/setprec.o setprec.c

$(OBJ_DIR)/sin.o: sin.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sin.o sin.c

$(OBJ_DIR)/sincos.o: sincos.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sincos.o sincos.c

$(OBJ_DIR)/sindg.o: sindg.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sindg.o sindg.c

$(OBJ_DIR)/sinh.o: sinh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sinh.o sinh.c

$(OBJ_DIR)/sqrt.o: sqrt.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sqrt.o sqrt.c

$(OBJ_DIR)/tan.o: tan.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/tan.o tan.c

$(OBJ_DIR)/tandg.o: tandg.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/tandg.o tandg.c

$(OBJ_DIR)/tanh.o: tanh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/tanh.o tanh.c

$(OBJ_DIR)/unity.o: unity.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/unity.o unity.c


clean:
	rm -f $(OBJS)

