
PROG = x2dhf
XHF = ../../tests/xhf 

all: h2 h2a h2b 
	$(MAKE) clean

h2: h2.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@-1 | tee $@.lst

h2a: h2a.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@-1 | tee $@.lst

h2b: h2b.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@-1 | tee $@.lst

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

