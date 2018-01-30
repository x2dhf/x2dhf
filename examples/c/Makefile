
PROG = x2dhf
XHF = ../../tests/xhf 

all: c c1d2 c3p1
	$(MAKE) clean

c: c.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@-1 | tee $@.lst

c1d2: c1d2.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@-1 | tee $@.lst

c3p1: c3p1.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@-1 | tee $@.lst


clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#


#

