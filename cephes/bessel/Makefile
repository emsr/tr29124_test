
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

INC_DIR = ../cmath
OBJ_DIR = obj

OBJS = \
  $(OBJ_DIR)/airy.o \
  $(OBJ_DIR)/hyp2f1.o \
  $(OBJ_DIR)/hyperg.o \
  $(OBJ_DIR)/i0.o \
  $(OBJ_DIR)/i1.o \
  $(OBJ_DIR)/iv.o \
  $(OBJ_DIR)/j0.o \
  $(OBJ_DIR)/j1.o \
  $(OBJ_DIR)/jn.o \
  $(OBJ_DIR)/jv.o \
  $(OBJ_DIR)/k0.o \
  $(OBJ_DIR)/k1.o \
  $(OBJ_DIR)/kn.o \
  $(OBJ_DIR)/psi.o \
  $(OBJ_DIR)/struve.o \
  $(OBJ_DIR)/yn.o

libcephes_bessel.so: $(OBJS)
	$(CC) -shared -o libcephes_bessel.so $(OBJS) -L.. -lcephes_cmath -lcephes_cprob
	cp libcephes_bessel.so libcephes_bessel.dll
	cp libcephes_bessel.so libcephes_bessel.dll ..

$(OBJ_DIR)/airy.o: airy.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/airy.o airy.c

$(OBJ_DIR)/hyp2f1.o: hyp2f1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/hyp2f1.o hyp2f1.c

$(OBJ_DIR)/hyperg.o: hyperg.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/hyperg.o hyperg.c

$(OBJ_DIR)/i0.o: i0.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/i0.o i0.c

$(OBJ_DIR)/i1.o: i1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/i1.o i1.c

$(OBJ_DIR)/iv.o: iv.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/iv.o iv.c

$(OBJ_DIR)/j0.o: j0.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/j0.o j0.c

$(OBJ_DIR)/j1.o: j1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/j1.o j1.c

$(OBJ_DIR)/jn.o: jn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/jn.o jn.c

$(OBJ_DIR)/jv.o: jv.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/jv.o jv.c

$(OBJ_DIR)/k0.o: k0.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/k0.o k0.c

$(OBJ_DIR)/k1.o: k1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/k1.o k1.c

$(OBJ_DIR)/kn.o: kn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/kn.o kn.c

$(OBJ_DIR)/psi.o: psi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/psi.o psi.c

$(OBJ_DIR)/struve.o: struve.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/struve.o struve.c

$(OBJ_DIR)/yn.o: yn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/yn.o yn.c


clean:
	rm -f $(OBJS)

