
CC = $(HOME)/bin/bin/gcc -std=c11

OBJS = bdtr.o btdtr.o chdtr.o const.o drand.o expx2.o fdtr.o gamma.o gdtr.o igam.o igami.o incbet.o incbi.o kolmogorov.o mtherr.o nbdtr.o ndtr.o ndtri.o pdtr.o polevl.o stdtr.o unity.o

INC_DIR = .

libcephes_cprob.so: $(OBJS)
	$(CC) -shared -o libcephes_cprob.so $(OBJS) -lm

bdtr.o: bdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o bdtr.o bdtr.c

btdtr.o: btdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o btdtr.o btdtr.c

chdtr.o: chdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o chdtr.o chdtr.c

const.o: const.c
	$(CC) -c -fPIC -I$(INC_DIR) -o const.o const.c

drand.o: drand.c
	$(CC) -c -fPIC -I$(INC_DIR) -o drand.o drand.c

expx2.o: expx2.c
	$(CC) -c -fPIC -I$(INC_DIR) -o expx2.o expx2.c

fdtr.o: fdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o fdtr.o fdtr.c

gamma.o: gamma.c
	$(CC) -c -fPIC -I$(INC_DIR) -o gamma.o gamma.c

gdtr.o: gdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o gdtr.o gdtr.c

igam.o: igam.c
	$(CC) -c -fPIC -I$(INC_DIR) -o igam.o igam.c

igami.o: igami.c
	$(CC) -c -fPIC -I$(INC_DIR) -o igami.o igami.c

incbet.o: incbet.c
	$(CC) -c -fPIC -I$(INC_DIR) -o incbet.o incbet.c

incbi.o: incbi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o incbi.o incbi.c

kolmogorov.o: kolmogorov.c
	$(CC) -c -fPIC -I$(INC_DIR) -o kolmogorov.o kolmogorov.c

mtherr.o: mtherr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o mtherr.o mtherr.c

nbdtr.o: nbdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o nbdtr.o nbdtr.c

ndtr.o: ndtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o ndtr.o ndtr.c

ndtri.o: ndtri.c
	$(CC) -c -fPIC -I$(INC_DIR) -o ndtri.o ndtri.c

pdtr.o: pdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o pdtr.o pdtr.c

polevl.o: polevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o polevl.o polevl.c

stdtr.o: stdtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o stdtr.o stdtr.c

unity.o: unity.c
	$(CC) -c -fPIC -I$(INC_DIR) -o unity.o unity.c


clean:
	rm -f *.o
