%.o : %.f90
	$(F77) $(F77OPTS) -c $*.f90

all: Res

SRCLIBS= 

DRIVERS= driver
RESULTS= Res

Objs= driver.o air.o cons.o parab.o
driver: $(Objs)
	$(F77) $(F77OPTS) -o driver $(Objs) $(SRCLIBS)
Res: driver
	./driver >Res

air.o: cons.o

cons.o:

driver.o: cons.o parab.o

parab.o: cons.o air.o

diffres: Res res
	echo "Differences in results from driver"
	$(DIFF)  Res res

clean: 
	rm -rf *.o $(DRIVERS) $(CLEANUP) $(RESULTS)
