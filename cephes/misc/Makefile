
CC = $(HOME)/bin/bin/gcc -std=c11

INC_DIR = .

SRCS = \
  beta.c \
  chbevl.c \
  dawsn.c \
  ei.c \
  expn.c \
  fac.c \
  fresnl.c \
  mtherr.c \
  planck.c \
  polevl.c \
  polylog.c \
  psi.c \
  revers.c \
  rgamma.c \
  shichi.c \
  sici.c \
  simpsn.c \
  spence.c \
  zeta.c \
  zetac.c

OBJS = \
  beta.o \
  chbevl.o \
  dawsn.o \
  ei.o \
  expn.o \
  fac.o \
  fresnl.o \
  mtherr.o \
  planck.o \
  polevl.o \
  polylog.o \
  psi.o \
  revers.o \
  rgamma.o \
  shichi.o \
  sici.o \
  simpsn.o \
  spence.o \
  zeta.o \
  zetac.o

libcephes_misc.so: cephes_misc.h $(OBJS)
	$(CC) -shared -o libcephes_misc.so $(OBJS) -lm

beta.o: beta.c
	$(CC) -c -fPIC -I$(INC_DIR) -o beta.o beta.c

chbevl.o: chbevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o chbevl.o chbevl.c

dawsn.o: dawsn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o dawsn.o dawsn.c

ei.o: ei.c
	$(CC) -c -fPIC -I$(INC_DIR) -o ei.o ei.c

expn.o: expn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o expn.o expn.c

fac.o: fac.c
	$(CC) -c -fPIC -I$(INC_DIR) -o fac.o fac.c

fresnl.o: fresnl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o fresnl.o fresnl.c

mtherr.o: mtherr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o mtherr.o mtherr.c

planck.o: planck.c
	$(CC) -c -fPIC -I$(INC_DIR) -o planck.o planck.c

polevl.o: polevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o polevl.o polevl.c

polylog.o: polylog.c
	$(CC) -c -fPIC -I$(INC_DIR) -o polylog.o polylog.c

psi.o: psi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o psi.o psi.c

revers.o: revers.c
	$(CC) -c -fPIC -I$(INC_DIR) -o revers.o revers.c

rgamma.o: rgamma.c
	$(CC) -c -fPIC -I$(INC_DIR) -o rgamma.o rgamma.c

shichi.o: shichi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o shichi.o shichi.c

sici.o: sici.c
	$(CC) -c -fPIC -I$(INC_DIR) -o sici.o sici.c

simpsn.o: simpsn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o simpsn.o simpsn.c

spence.o: spence.c
	$(CC) -c -fPIC -I$(INC_DIR) -o spence.o spence.c

zeta.o: zeta.c
	$(CC) -c -fPIC -I$(INC_DIR) -o zeta.o zeta.c

zetac.o: zetac.c
	$(CC) -c -fPIC -I$(INC_DIR) -o zetac.o zetac.c

