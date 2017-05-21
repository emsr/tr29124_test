OBJECTS = nrcheck.o nrlinwww.o nrutil.o nrframe.o
SOURCE  = nrcheck.c nrdemo.c nrutil.c nrframe.c nrlinwww.c

LIBS       =  -lm

CC = cc 

nrdemo		:	nrdemo.o $(OBJECTS)
	$(CC) -o nrdemo nrdemo.o $(OBJECTS) $(LIBS)

clean		:	
	rm *.o
