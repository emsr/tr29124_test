
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

INC_DIR = .
OBJ_DIR = obj

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
  $(OBJ_DIR)/beta.o \
  $(OBJ_DIR)/chbevl.o \
  $(OBJ_DIR)/dawsn.o \
  $(OBJ_DIR)/ei.o \
  $(OBJ_DIR)/expn.o \
  $(OBJ_DIR)/fac.o \
  $(OBJ_DIR)/fresnl.o \
  $(OBJ_DIR)/mtherr.o \
  $(OBJ_DIR)/planck.o \
  $(OBJ_DIR)/polevl.o \
  $(OBJ_DIR)/polylog.o \
  $(OBJ_DIR)/psi.o \
  $(OBJ_DIR)/revers.o \
  $(OBJ_DIR)/rgamma.o \
  $(OBJ_DIR)/shichi.o \
  $(OBJ_DIR)/sici.o \
  $(OBJ_DIR)/simpsn.o \
  $(OBJ_DIR)/spence.o \
  $(OBJ_DIR)/zeta.o \
  $(OBJ_DIR)/zetac.o

libcephes_misc.so: cephes_misc.h $(OBJS)
	$(CC) -shared -o libcephes_misc.so $(OBJS) -L.. -lcephes_cmath -lcephes_cprob -lcephes_polyn
	cp libcephes_misc.so libcephes_misc.dll
	cp libcephes_misc.so libcephes_misc.dll ..

$(OBJ_DIR)/beta.o: beta.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/beta.o beta.c

$(OBJ_DIR)/chbevl.o: chbevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/chbevl.o chbevl.c

$(OBJ_DIR)/dawsn.o: dawsn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/dawsn.o dawsn.c

$(OBJ_DIR)/ei.o: ei.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ei.o ei.c

$(OBJ_DIR)/expn.o: expn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/expn.o expn.c

$(OBJ_DIR)/fac.o: fac.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/fac.o fac.c

$(OBJ_DIR)/fresnl.o: fresnl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/fresnl.o fresnl.c

$(OBJ_DIR)/mtherr.o: mtherr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mtherr.o mtherr.c

$(OBJ_DIR)/planck.o: planck.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/planck.o planck.c

$(OBJ_DIR)/polevl.o: polevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polevl.o polevl.c

$(OBJ_DIR)/polylog.o: polylog.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polylog.o polylog.c

$(OBJ_DIR)/psi.o: psi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/psi.o psi.c

$(OBJ_DIR)/revers.o: revers.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/revers.o revers.c

$(OBJ_DIR)/rgamma.o: rgamma.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/rgamma.o rgamma.c

$(OBJ_DIR)/shichi.o: shichi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/shichi.o shichi.c

$(OBJ_DIR)/sici.o: sici.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sici.o sici.c

$(OBJ_DIR)/simpsn.o: simpsn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/simpsn.o simpsn.c

$(OBJ_DIR)/spence.o: spence.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/spence.o spence.c

$(OBJ_DIR)/zeta.o: zeta.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/zeta.o zeta.c

$(OBJ_DIR)/zetac.o: zetac.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/zetac.o zetac.c


clean:
	rm -f $(OBJS)

