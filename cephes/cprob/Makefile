
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

# Let cmath have const.
#  $(OBJ_DIR)/const.o

INC_DIR = .
OBJ_DIR = obj

OBJS = \
  $(OBJ_DIR)/bdtr.o \
  $(OBJ_DIR)/btdtr.o \
  $(OBJ_DIR)/chdtr.o \
  $(OBJ_DIR)/drand.o \
  $(OBJ_DIR)/expx2.o \
  $(OBJ_DIR)/fdtr.o \
  $(OBJ_DIR)/gamma.o \
  $(OBJ_DIR)/gdtr.o \
  $(OBJ_DIR)/igam.o \
  $(OBJ_DIR)/igami.o \
  $(OBJ_DIR)/incbet.o \
  $(OBJ_DIR)/incbi.o \
  $(OBJ_DIR)/kolmogorov.o \
  $(OBJ_DIR)/mtherr.o \
  $(OBJ_DIR)/nbdtr.o \
  $(OBJ_DIR)/ndtr.o \
  $(OBJ_DIR)/ndtri.o \
  $(OBJ_DIR)/pdtr.o \
  $(OBJ_DIR)/polevl.o \
  $(OBJ_DIR)/stdtr.o \
  $(OBJ_DIR)/unity.o

INC_DIR = .
OBJ_DIR = obj

libcephes_cprob.so: $(OBJS)
	$(CC) -shared -o libcephes_cprob.so $(OBJS) -L.. -lcephes_cmath
	cp libcephes_cprob.so libcephes_cprob.dll
	cp libcephes_cprob.so libcephes_cprob.dll ..

$(OBJ_DIR)/bdtr.o: bdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/bdtr.o bdtr.c

$(OBJ_DIR)/btdtr.o: btdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/btdtr.o btdtr.c

$(OBJ_DIR)/chdtr.o: chdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/chdtr.o chdtr.c

$(OBJ_DIR)/const.o: const.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/const.o const.c

$(OBJ_DIR)/drand.o: drand.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/drand.o drand.c

$(OBJ_DIR)/expx2.o: expx2.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/expx2.o expx2.c

$(OBJ_DIR)/fdtr.o: fdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/fdtr.o fdtr.c

$(OBJ_DIR)/gamma.o: gamma.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/gamma.o gamma.c

$(OBJ_DIR)/gdtr.o: gdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/gdtr.o gdtr.c

$(OBJ_DIR)/igam.o: igam.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/igam.o igam.c

$(OBJ_DIR)/igami.o: igami.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/igami.o igami.c

$(OBJ_DIR)/incbet.o: incbet.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/incbet.o incbet.c

$(OBJ_DIR)/incbi.o: incbi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/incbi.o incbi.c

$(OBJ_DIR)/kolmogorov.o: kolmogorov.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/kolmogorov.o kolmogorov.c

$(OBJ_DIR)/mtherr.o: mtherr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mtherr.o mtherr.c

$(OBJ_DIR)/nbdtr.o: nbdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/nbdtr.o nbdtr.c

$(OBJ_DIR)/ndtr.o: ndtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ndtr.o ndtr.c

$(OBJ_DIR)/ndtri.o: ndtri.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ndtri.o ndtri.c

$(OBJ_DIR)/pdtr.o: pdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/pdtr.o pdtr.c

$(OBJ_DIR)/polevl.o: polevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polevl.o polevl.c

$(OBJ_DIR)/stdtr.o: stdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/stdtr.o stdtr.c

$(OBJ_DIR)/unity.o: unity.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/unity.o unity.c


clean:
	rm -f $(OBJS)
