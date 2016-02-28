
CC = $(HOME)/bin/bin/gcc

OBJS = airy.o hyp2f1.o hyperg.o i0.o i1.o iv.o j0.o j1.o jn.o jv.o k0.o k1.o kn.o psi.o struve.o yn.o

INC_DIR = ../128bit

libcephes_bessel.so: $(OBJS)
	$(CC) -shared -o libcephes_bessel.so $(OBJS)

airy.o: airy.c
	$(CC) -c -fPIC -I$(INC_DIR) -o airy.o airy.c

hyp2f1.o: hyp2f1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o hyp2f1.o hyp2f1.c

hyperg.o: hyperg.c
	$(CC) -c -fPIC -I$(INC_DIR) -o hyperg.o hyperg.c

i0.o: i0.c
	$(CC) -c -fPIC -I$(INC_DIR) -o i0.o i0.c

i1.o: i1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o i1.o i1.c

iv.o: iv.c
	$(CC) -c -fPIC -I$(INC_DIR) -o iv.o iv.c

j0.o: j0.c
	$(CC) -c -fPIC -I$(INC_DIR) -o j0.o j0.c

j1.o: j1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o j1.o j1.c

jn.o: jn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o jn.o jn.c

jv.o: jv.c
	$(CC) -c -fPIC -I$(INC_DIR) -o jv.o jv.c

k0.o: k0.c
	$(CC) -c -fPIC -I$(INC_DIR) -o k0.o k0.c

k1.o: k1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o k1.o k1.c

kn.o: kn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o kn.o kn.c

psi.o: psi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o psi.o psi.c

struve.o: struve.c
	$(CC) -c -fPIC -I$(INC_DIR) -o struve.o struve.c

yn.o: yn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o yn.o yn.c


clean:
	rm -f *.o

