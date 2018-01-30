
PROG = x2dhf
XHF = ../../tests/xhf 

all: be-1 be-2
	$(MAKE) clean

be-1: be-1.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ be-b88  | tee $@.lst

be-2: be-2.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ be-b88  be-b88+lyp | tee $@.lst

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

