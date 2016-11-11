
CC = $(HOME)/bin/bin/gcc

OBJS = euclid.o \
 polmisc.o \
 polrt.o \
 polyn.o \
 polyr.o \
 revers.o

INC_DIR = ../128bit

libcephes_polyn.so: $(OBJS)
	$(CC) -shared -o libcephes_polyn.so $(OBJS)
	cp libcephes_polyn.so libcephes_polyn.dll
	cp libcephes_polyn.so libcephes_polyn.dll ..

euclid.o: euclid.c
	$(CC) -c -fPIC -I$(INC_DIR) -o euclid.o euclid.c

polmisc.o: polmisc.c
	$(CC) -c -fPIC -I$(INC_DIR) -o polmisc.o polmisc.c

polrt.o: polrt.c
	$(CC) -c -fPIC -I$(INC_DIR) -o polrt.o polrt.c

polyn.o: polyn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o polyn.o polyn.c

polyr.o: polyr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o polyr.o polyr.c

revers.o: revers.c
	$(CC) -c -fPIC -I$(INC_DIR) -o revers.o revers.c

