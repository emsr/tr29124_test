include makefile.inc

OBJS = driver.o bsslr.o mathur.o mbsslr.o mcnr.o mmathur.o

all: res1 res2

res1: $(OBJS) data1
	$(CPP) $(CPPOPTS) -o driver $(OBJS)
	driver < data1 > res1

res2: $(OBJS) data2
	$(CPP) $(CPPOPTS) -o driver $(OBJS)
	driver < data2 > res2

driver.o: driver.cpp $(INCSRC)/bsslr.cpp $(INCSRC)/mathur.cpp \
          $(INCSRC)/mbsslr.cpp $(INCSRC)/mcnr.cpp $(INCSRC)/mmathur.cpp
	$(CPP) $(CPPOPTS) -c -I$(INCSRC) driver.cpp

bsslr.o: $(INCSRC)/bsslr.cpp
	$(CPP) $(CPPOPTS) -c $(INCSRC)/bsslr.cpp

mathur.o: $(INCSRC)/mathur.cpp
	$(CPP) $(CPPOPTS) -c $(INCSRC)/mathur.cpp

mbsslr.o: $(INCSRC)/mbsslr.cpp
	$(CPP) $(CPPOPTS) -c $(INCSRC)/mbsslr.cpp

mcnr.o: $(INCSRC)/mcnr.cpp
	$(CPP) $(CPPOPTS) -c $(INCSRC)/mcnr.cpp

mmathur.o: $(INCSRC)/mmathur.cpp
	$(CPP) $(CPPOPTS) -c $(INCSRC)/mmathur.cpp
