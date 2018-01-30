
PROG = x2dhf
XHF = ../../tests/xhf 

all: bf_g94 bf_g94l

bf_g94: bf_g94.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || (cp $@.pun gaussian.pun; cp $@.out gaussian.out; $(XHF) $@ $@-1 | tee $@.lst)

bf_g94l: bf_g94l.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || (cp $@.out gaussianc.out; $(XHF) $@ $@-1 | tee $@.lst)
	$(MAKE) clean

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) gaussian.* gaussianc.out 2> /dev/null

stop:
	$(XHF) stop

#
