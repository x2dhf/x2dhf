
PROG = x2dhf
XHF = ../../tests/xhf 

all: c2a c2b c2c
	$(MAKE) clean

c2a: c2a.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || cp c2.pun gaussian.pun; cp c2.out gaussian.out; $(XHF) $@ $@-1 | tee $@.lst

c2b: c2b.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ c2a-1 $@-1 | tee $@.lst

c2c: c2c.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ c2b-1 $@-1 | tee $@.lst


clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) gaussian.* gaussianc.out 2> /dev/null

stop:
	$(XHF) stop

#

