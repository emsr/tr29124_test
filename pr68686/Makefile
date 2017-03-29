tgammaq-bug: tgammaq-bug.c
	$(HOME)/bin/bin/gcc -g -o tgammaq-bug tgammaq-bug.c -lquadmath

pr68686: pr68686.c
	$(HOME)/bin/bin/gcc -g -o pr68686 pr68686.c -lquadmath

test:
	PATH=$(HOME)/bin/bin:$$PATH ./tgammaq-bug

check:
	PATH=$(HOME)/bin/bin:$$PATH ./pr68686
