all: tgammaq-bug tgammaq-test pr68686

tgammaq-bug: tgammaq-bug.c
	$(HOME)/bin/bin/gcc -g -o tgammaq-bug tgammaq-bug.c -lquadmath

tgammaq-test: tgammaq-test.c
	$(HOME)/bin/bin/gcc -g -o tgammaq-test tgammaq-test.c -lquadmath -lm

pr68686: pr68686.c
	$(HOME)/bin/bin/gcc -g -o pr68686 pr68686.c -lquadmath

test:
	./tgammaq-bug > tgammaq-bug.buggy
	LD_LIBRARY_PATH=$(HOME)/bin/lib64:$(LD_LIBRARY_PATH) ./tgammaq-bug > tgammaq-bug.patched
	rm -f tgamma.dat
	rm -f tgammaq.dat
	
	./tgammaq-test
	mv tgamma.dat tgamma-buggy.dat
	mv tgammaq.dat tgammaq-buggy.dat
	LD_LIBRARY_PATH=$(HOME)/bin/lib64:$(LD_LIBRARY_PATH) ./tgammaq-test
	mv tgamma.dat tgamma-patched.dat
	mv tgammaq.dat tgammaq-patched.dat

check:
	LD_LIBRARY_PATH=$(HOME)/bin/lib64:$(LD_LIBRARY_PATH) ./pr68686 2> pr68686-patched.txt
	./pr68686 2> pr68686-buggy.txt
