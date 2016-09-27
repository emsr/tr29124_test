FFLAGS= -sloppy -C -d5 -g -temp=/tmp -u
FC=epcf90

OBJS = driver.o src.o airy.f port.f

res : $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) 
	a.out
