
CC = $(HOME)/bin/bin/gcc -std=c11

INC_DIR = .

SRCS = \
  cmplx.c \
  const.c \
  ellf.c \
  ellie.c \
  ellik.c \
  ellpe.c \
  ellpj.c \
  ellpk.c \
  mtherr.c \
  polevl.c

# This is a main program.
#  ellf.o \

OBJS = \
  cmplx.o \
  const.o \
  ellie.o \
  ellik.o \
  ellpe.o \
  ellpj.o \
  ellpk.o \
  mtherr.o \
  polevl.o

libcephes_ellf.so: cephes_ellf.h $(OBJS)
	$(CC) -shared -o libcephes_ellf.so $(OBJS) -lm
	cp libcephes_ellf.so libcephes_ellf.dll
	cp libcephes_ellf.so libcephes_ellf.dll ..

cmplx.o: cmplx.c
	$(CC) -c -fPIC -I$(INC_DIR) -o cmplx.o cmplx.c

const.o: const.c
	$(CC) -c -fPIC -I$(INC_DIR) -o const.o const.c

ellf.o: ellf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o ellf.o ellf.c

ellie.o: ellie.c
	$(CC) -c -fPIC -I$(INC_DIR) -o ellie.o ellie.c

ellik.o: ellik.c
	$(CC) -c -fPIC -I$(INC_DIR) -o ellik.o ellik.c

ellpe.o: ellpe.c
	$(CC) -c -fPIC -I$(INC_DIR) -o ellpe.o ellpe.c

ellpj.o: ellpj.c
	$(CC) -c -fPIC -I$(INC_DIR) -o ellpj.o ellpj.c

ellpk.o: ellpk.c
	$(CC) -c -fPIC -I$(INC_DIR) -o ellpk.o ellpk.c

mtherr.o: mtherr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o mtherr.o mtherr.c

polevl.o: polevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o polevl.o polevl.c

