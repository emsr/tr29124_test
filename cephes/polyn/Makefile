
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

INC_DIR = ../cmath
OBJ_DIR = obj

# Rational polynomials.
#  $(OBJ_DIR)/polyr.o

OBJS = \
  $(OBJ_DIR)/euclid.o \
  $(OBJ_DIR)/polmisc.o \
  $(OBJ_DIR)/polrt.o \
  $(OBJ_DIR)/polyn.o \
  $(OBJ_DIR)/revers.o

libcephes_polyn.so: $(OBJS)
	$(CC) -shared -o libcephes_polyn.so $(OBJS) -L.. -lcephes_cmath
	cp libcephes_polyn.so libcephes_polyn.dll
	cp libcephes_polyn.so libcephes_polyn.dll ..

$(OBJ_DIR)/euclid.o: euclid.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/euclid.o euclid.c

$(OBJ_DIR)/polmisc.o: polmisc.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polmisc.o polmisc.c

$(OBJ_DIR)/polrt.o: polrt.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polrt.o polrt.c

$(OBJ_DIR)/polyn.o: polyn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polyn.o polyn.c

$(OBJ_DIR)/polyr.o: polyr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polyr.o polyr.c

$(OBJ_DIR)/revers.o: revers.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/revers.o revers.c


clean:
	rm -f $(OBJS)
