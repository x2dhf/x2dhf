
PROG = x2dhf
XHF = ../../tests/xhf 

all: fh-0 fh-m1 fh-m2 fh-p1 fh-p2 fh-elprop 
	$(MAKE) clean

fh-0: fh-0.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@-1 | tee $@.lst

fh-m1: fh-m1.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ fh-0-1 $@-1 | tee $@.lst

fh-m2: fh-m2.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ fh-m1-1 $@-1 | tee $@.lst

fh-p1: fh-p1.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ fh-0-1  $@-1 | tee $@.lst

fh-p2: fh-p2.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ fh-p1-1 $@-1 | tee $@.lst

fh-elprop: 
	./fh-elprop.sh

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

