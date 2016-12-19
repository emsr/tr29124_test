
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

INC_DIR = .
OBJ_DIR = obj

SRCS = \
  cgamma.c \
  cmplx.c \
  cpolylog.c

OBJS = \
  $(OBJ_DIR)/cmplx.o \
  $(OBJ_DIR)/cgamma.o \
  $(OBJ_DIR)/cpolylog.o

libcephes_cmplx.so: $(OBJS)
	$(CC) -shared -o libcephes_cmplx.so $(OBJS) -L.. -lcephes_cmath -lcephes_misc
	cp libcephes_cmplx.so libcephes_cmplx.dll
	cp libcephes_cmplx.so libcephes_cmplx.dll ..

$(OBJ_DIR)/cmplx.o: cmplx.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cmplx.o cmplx.c

$(OBJ_DIR)/cgamma.o: cgamma.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cgamma.o cgamma.c

$(OBJ_DIR)/cpolylog.o: cpolylog.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cpolylog.o cpolylog.c
