LD_FLAGS= -lm
CPPFLAGS= -I../Src/

all: driver test_boundary Results.dat

driver: driver.o ../Src/wright.o
	$(CC) -o $@ $^ $(LD_FLAGS)

test_boundary: test_boundary.o ../Src/wright.o
	$(CC) -o $@ $^ $(LD_FLAGS)

Results.dat: test_boundary
	./test_boundary

diff_boundary: test_boundary
	echo "Differences in results from driver:"
	diff Results.dat results.dat

.PHONY:clean
clean:
	rm -f *.o driver test_boundary Results.dat
