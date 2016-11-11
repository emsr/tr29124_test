
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

INC_DIR = .
OBJ_DIR = obj

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
#  ellf.o
# Let cmath have const.
#  $(OBJ_DIR)/const.o

OBJS = \
  $(OBJ_DIR)/cmplx.o \
  $(OBJ_DIR)/ellie.o \
  $(OBJ_DIR)/ellik.o \
  $(OBJ_DIR)/ellpe.o \
  $(OBJ_DIR)/ellpj.o \
  $(OBJ_DIR)/ellpk.o \
  $(OBJ_DIR)/mtherr.o \
  $(OBJ_DIR)/polevl.o

libcephes_ellf.so: cephes_ellf.h $(OBJS)
	$(CC) -shared -o libcephes_ellf.so $(OBJS) -L.. -lcephes_cmath
	cp libcephes_ellf.so libcephes_ellf.dll
	cp libcephes_ellf.so libcephes_ellf.dll ..

$(OBJ_DIR)/cmplx.o: cmplx.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cmplx.o cmplx.c

$(OBJ_DIR)/const.o: const.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/const.o const.c

$(OBJ_DIR)/ellf.o: ellf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellf.o ellf.c

$(OBJ_DIR)/ellie.o: ellie.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellie.o ellie.c

$(OBJ_DIR)/ellik.o: ellik.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellik.o ellik.c

$(OBJ_DIR)/ellpe.o: ellpe.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellpe.o ellpe.c

$(OBJ_DIR)/ellpj.o: ellpj.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellpj.o ellpj.c

$(OBJ_DIR)/ellpk.o: ellpk.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ellpk.o ellpk.c

$(OBJ_DIR)/mtherr.o: mtherr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mtherr.o mtherr.c

$(OBJ_DIR)/polevl.o: polevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polevl.o polevl.c


clean:
	rm -f $(OBJS)

